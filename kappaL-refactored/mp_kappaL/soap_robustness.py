"""Step 10：结构表征稳健性 —— 晶胞标准化 sanity check + 几何/组成权重敏感性。

1. 晶胞标准化：对 100 个随机材料分别用原始胞与 SpacegroupAnalyzer primitive standard
   胞算 SOAP，报告余弦相似度分布（标准化后应接近 1）。
2. 权重敏感性：Structure = w*几何 + (1-w)*组成，w ∈ {0,0.25,0.5,0.75,1}，
   跑 Structure↔kL 重叠 + 距离 Spearman，看结论随 w 的变化。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import config
from mp_kappaL.data_utils import load_aligned
from mp_kappaL import metrics
from mp_kappaL.build_views import parse_structure_ase, parse_structure_pymatgen, standardize
from graph_utils import hellinger_distance, soap_distance
from dscribe.descriptors import SOAP


def _soap_of(atoms_list):
    soap = SOAP(species=config.SOAP_SPECIES, periodic=True, r_cut=config.SOAP_R_CUT,
                n_max=config.SOAP_N_MAX, l_max=config.SOAP_L_MAX,
                sigma=config.SOAP_SIGMA, average=config.SOAP_AVERAGE)
    return np.asarray(soap.create(atoms_list, n_jobs=-1), dtype=np.float32)


def standardization_sanity(n=100):
    """100 个材料：原始胞 vs primitive standard 胞的 SOAP 余弦相似度。"""
    el = json.load(open(config.RAW_DIR / "elasticity_all.json"))
    rng = np.random.RandomState(config.SEED)
    cands = [x for x in el if x.get("structure")]
    pick = rng.choice(cands, size=n, replace=False)
    cos_random, cos_primitive, cos_conventional = [], [], []
    n_ok = 0
    for x in pick:
        st = x["structure"]
        try:
            prim = standardize(st)
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            conv = SpacegroupAnalyzer(
                parse_structure_pymatgen(st), symprec=0.1).get_conventional_standard_structure()
            from ase import Atoms
            a_orig = parse_structure_ase(st)
            a_prim = Atoms(symbols="X" * len(prim), positions=prim.cart_coords,
                           cell=prim.lattice.matrix, pbc=True)
            a_conv = Atoms(symbols="X" * len(conv), positions=conv.cart_coords,
                           cell=conv.lattice.matrix, pbc=True)
            v_orig, v_prim, v_conv = _soap_of([a_orig, a_prim, a_conv])
            v_orig = v_orig / (np.linalg.norm(v_orig) + 1e-12)
            v_prim = v_prim / (np.linalg.norm(v_prim) + 1e-12)
            v_conv = v_conv / (np.linalg.norm(v_conv) + 1e-12)
            cos_primitive.append(float(v_orig @ v_prim))
            cos_conventional.append(float(v_orig @ v_conv))
            # 也计算两个不同原始胞之间的基线相似度（随机向量）
            w = rng.standard_normal(v_orig.shape[0]); w /= np.linalg.norm(w)
            cos_random.append(float(v_orig @ w))
            n_ok += 1
        except Exception as e:
            continue
    cos_primitive = np.array(cos_primitive)
    cos_conventional = np.array(cos_conventional)
    cos_random = np.array(cos_random)
    summ = pd.DataFrame({
        "metric": ["orig_vs_primitive_cosine", "orig_vs_conventional_cosine",
                   "orig_vs_random_cosine"],
        "mean": [cos_primitive.mean(), cos_conventional.mean(), cos_random.mean()],
        "min": [cos_primitive.min(), cos_conventional.min(), cos_random.min()],
        "p50": [np.median(cos_primitive), np.median(cos_conventional),
                np.median(cos_random)],
    })
    print(f"标准化 sanity check: {n_ok}/{n} 成功")
    print(summ.to_string(index=False))
    summ.to_csv(config.PROC_DIR / "soap_standardization_sanity.csv", index=False)
    return summ


def weight_sensitivity():
    """Structure = w*geo + (1-w)*comp 的 w 敏感性。"""
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    soap_geo = feats["soap_geo"].astype(np.float32)
    comp_frac = feats["comp_frac"].astype(np.float32)
    if len(meta) > config.WEIGHT_AUDIT_N:
        audit_rng = np.random.RandomState(config.SEED)
        take = np.sort(audit_rng.choice(len(meta), config.WEIGHT_AUDIT_N, replace=False))
        meta = meta.iloc[take].reset_index(drop=True)
        soap_geo = soap_geo[take]
        comp_frac = comp_frac[take]
    d_geo_raw = soap_distance(soap_geo).astype(np.float32)
    d_comp_raw = hellinger_distance(comp_frac).astype(np.float32)

    ids = meta["material_id"].to_numpy()
    rng = np.random.RandomState(config.SEED)

    rows = []
    for normalization in ["p95", "max"]:
        q = 0.95 if normalization == "p95" else 1.0
        geo_scale = np.quantile(d_geo_raw[d_geo_raw > 0], q)
        comp_scale = np.quantile(d_comp_raw[d_comp_raw > 0], q)
        d_geo = np.clip(d_geo_raw / geo_scale, 0, 1)
        d_comp = np.clip(d_comp_raw / comp_scale, 0, 1)
        for target in ["snyder_acoustic", "clarke"]:
            y = np.log10(meta[target].values.astype(float))
            d_kL = np.abs(y[:, None] - y[None, :]).astype(np.float32)
            nnK = metrics.knn_neighbor_matrix(d_kL, config.K, ids=ids)
            for w in config.W_GEO_GRID:
                d_struct = (w * d_geo + (1 - w) * d_comp).astype(np.float32)
                nnS = metrics.knn_neighbor_matrix(d_struct, config.K, ids=ids)
                r = metrics.crossview_overlap(nnS, nnK, config.K, rng, n_perm=100, n_boot=200)
                ip = rng.randint(0, len(ids), 50_000)
                jp = rng.randint(0, len(ids), 50_000)
                same = ip == jp
                while same.any():
                    jp[same] = rng.randint(0, len(ids), int(same.sum()))
                    same = ip == jp
                rho = stats.spearmanr(d_struct[ip, jp], d_kL[ip, jp]).statistic
                rows.append({"normalization": normalization, "target": target,
                             "w_geo": w, "overlap": r["overlap"], "z": r["z"],
                             "enrichment": r["enrichment"], "spearman": float(rho)})
                print(f"  norm={normalization} target={target} w_geo={w:.2f}: "
                      f"overlap={r['overlap']:.4f} Spearman={rho:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(config.PROC_DIR / "weight_sensitivity.csv", index=False)
    print("saved processed/weight_sensitivity.csv")
    return df


def main():
    standardization_sanity(n=100)
    print()
    print("=== 权重敏感性 ===")
    weight_sensitivity()


if __name__ == "__main__":
    main()
