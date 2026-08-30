"""MP 数据：构建 structure / elastic / kappa_L 视图（重构版，Step 2 / Step 10）。

变更（相对 v0）：
1. 物理合理性清洗在源头做一次（clean_records），输出 rejected_samples.csv 可审计。
2. 所有特征矩阵伴随显式 material_id 索引（row_index.npy）保存（Step 3）。
3. 默认先过 SpacegroupAnalyzer.get_primitive_standard_structure() 再做 SOAP；
   仅 ``--no-standardize`` 可显式关闭（Step 10）。
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config
from mp_kappaL.data_utils import clean_records

mp = config.MP_DIR
out = config.PROC_DIR
out.mkdir(parents=True, exist_ok=True)


def parse_structure_ase(st):
    """MP structure dict -> ASE Atoms（dummy X 物种，几何-only）。"""
    lat = np.array(st["lattice"]["matrix"], dtype=float)
    pos = np.array([s["xyz"] for s in st["sites"]], dtype=float)
    return Atoms(symbols="X" * len(pos), positions=pos, cell=lat, pbc=True)


def parse_structure_pymatgen(st):
    """MP structure dict -> pymatgen Structure（用于标准化）。"""
    from pymatgen.core import Structure
    return Structure.from_dict({
        "lattice": {"matrix": st["lattice"]["matrix"]},
        "sites": [{
            "species": [{"element": sp["element"], "occu": float(sp.get("occu", 1.0))}
                        for sp in s["species"]],
            "abc": s["abc"],
        } for s in st["sites"]],
    })


def standardize(st):
    """Step 10：SpacegroupAnalyzer 原始胞 -> primitive standard 胞。"""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    pmg = parse_structure_pymatgen(st)
    return SpacegroupAnalyzer(pmg, symprec=0.1).get_primitive_standard_structure()


def build_records(el):
    """从 elasticity_all.json 抽取有效样本（与 v0 相同的准入条件）。"""
    records = []
    for x in el:
        sv = x.get("sound_velocity") or {}
        tc = x.get("thermal_conductivity") or {}
        bm = x.get("bulk_modulus") or {}
        gm = x.get("shear_modulus") or {}
        st = x.get("structure")
        if (st and tc.get("clarke") and tc.get("clarke") > 0
                and bm.get("vrh") and gm.get("vrh")
                and x.get("debye_temperature") and x.get("density")):
            records.append({
                "material_id": x["material_id"],
                "formula": x.get("formula_pretty"),
                "structure": st,
                "clarke": float(tc["clarke"]),
                "cahill": float(tc.get("cahill", np.nan)),
                "snyder_total": float(sv.get("snyder_total", np.nan)),
                "snyder_acoustic": float(sv.get("snyder_acoustic", np.nan)),
                "bulk_vrh": float(bm["vrh"]),
                "shear_vrh": float(gm["vrh"]),
                "debye": float(x["debye_temperature"]),
                "density": float(x["density"]),
                "nsites": int(x.get("nsites") or 0),
                "v_long": float(sv.get("longitudinal", np.nan)),
                "v_trans": float(sv.get("transverse", np.nan)),
            })
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-standardize", dest="standardize", action="store_false",
                    help="关闭结构标准化（仅用于敏感性复验，不建议用于主流程）")
    ap.set_defaults(standardize=True)
    args = ap.parse_args()

    el = json.load(open(mp / "raw" / "elasticity_all.json"))
    records = build_records(el)
    print("准入 records（清洗前）:", len(records))

    meta = pd.DataFrame([{k: v for k, v in r.items() if k != "structure"} for r in records])
    meta_clean, rejected = clean_records(meta)
    n_before, n_after = len(meta), len(meta_clean)
    print(f"清洗: {n_before} -> {n_after}（剔除 {n_before - n_after} 条）")
    hits = rejected["rule"].value_counts() if len(rejected) else pd.Series(dtype=int)
    print("规则命中统计（独立计数，长表 rejected_samples.csv 一行=一次违规）:")
    print(hits.to_string())
    print("被剔除的唯一 material 数:", rejected["material_id"].nunique() if len(rejected) else 0)
    rejected.to_csv(out / "rejected_samples.csv", index=False)
    print("saved rejected_samples.csv:", len(rejected), "行（长表）")

    # 只保留清洗后样本的结构
    keep_ids = set(meta_clean["material_id"])
    clean_rec_list = [r for r in records if r["material_id"] in keep_ids]
    assert len(clean_rec_list) == len(meta_clean), "清洗后结构与 meta 不一致"

    # ---- 结构 -> ASE (dummy X) + 元素分数 ----
    elements_all = set()
    atoms_list = []
    comp_rows = []
    n_atom_changes = []
    for r in clean_rec_list:
        st = r["structure"]
        if args.standardize:
            prim = standardize(st)
            atoms_list.append(Atoms(symbols="X" * len(prim), positions=prim.cart_coords,
                                    cell=prim.lattice.matrix, pbc=True))
            n_atom_changes.append(len(prim) - len(st["sites"]))
        else:
            atoms_list.append(parse_structure_ase(st))
        frac = Counter()
        for s in st["sites"]:
            for sp in s["species"]:
                frac[sp["element"]] += float(sp.get("occu", 1.0))
        total = sum(frac.values()) or 1.0
        frac = {k: v / total for k, v in frac.items()}
        comp_rows.append(frac)
        elements_all.update(frac.keys())
    if args.standardize:
        print(f"标准化后原子数变化: 均 {np.mean(n_atom_changes):.2f}, "
              f"变多 {sum(v > 0 for v in n_atom_changes)}, "
              f"变少 {sum(v < 0 for v in n_atom_changes)}, "
              f"不变 {sum(v == 0 for v in n_atom_changes)}")

    elems = sorted(elements_all)
    F = np.zeros((len(clean_rec_list), len(elems)))
    for i, fr in enumerate(comp_rows):
        for e, v in fr.items():
            F[i, elems.index(e)] = v

    # ---- SOAP (dummy X) ----
    from dscribe.descriptors import SOAP
    soap = SOAP(species=config.SOAP_SPECIES, periodic=True, r_cut=config.SOAP_R_CUT,
                n_max=config.SOAP_N_MAX, l_max=config.SOAP_L_MAX,
                sigma=config.SOAP_SIGMA, average=config.SOAP_AVERAGE)
    geo = np.asarray(soap.create(atoms_list, n_jobs=-1), dtype=np.float32)
    print("SOAP geo shape:", geo.shape)

    # ---- 保存（特征矩阵伴随 row_index.npy）----
    row_ids = meta_clean["material_id"].to_numpy()
    np.save(out / "row_index.npy", row_ids)
    np.save(out / "soap_geo.npy", geo)
    np.save(out / "comp_frac.npy", F.astype(np.float32))
    np.save(out / "elem_basis.npy", np.array(elems))
    meta_clean.to_parquet(out / "views_meta.parquet", index=False)
    meta_clean.to_json(out / "views_meta.json")
    print(f"saved. clarke n={meta_clean['clarke'].notna().sum()} "
          f"snyder n={meta_clean['snyder_acoustic'].notna().sum()}")
    print("clarke range:", round(meta_clean["clarke"].min(), 4), "-",
          round(meta_clean["clarke"].max(), 4))


if __name__ == "__main__":
    main()
