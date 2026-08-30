"""MP 数据：kL 与单个物理描述符的直接相关（解释层）+ 图（重构版，Step 3）。

Step 3 修复：用 load_aligned 按 material_id 严格对齐特征矩阵与 meta，
删除了 v0 里「先 bulk_vrh<1000 过滤 meta 再用全长特征矩阵按行号索引」的错位 bug。
"""
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
from graph_utils import hellinger_distance, soap_distance

mp = config.MP_DIR


def main():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    soap_geo = feats["soap_geo"]
    comp_frac = feats["comp_frac"]
    elems = np.load(config.PROC_DIR / "elem_basis.npy", allow_pickle=True)

    from pymatgen.core import Element
    mass = np.array([float(Element(e).atomic_mass.real) for e in elems])
    avg_mass = (comp_frac * mass).sum(axis=1)
    meta["avg_mass"] = avg_mass
    meta["v_s_proxy"] = np.sqrt(meta["bulk_vrh"] * 1e9 / (meta["density"] * 1000))  # m/s

    logk = np.log10(meta["clarke"].values.astype(float))
    feats_list = ["debye", "bulk_vrh", "shear_vrh", "density", "avg_mass",
                  "v_long", "v_trans", "v_s_proxy"]

    print("=== Spearman correlation of log10(clarke) with descriptors ===")
    rows = []
    for c in feats_list:
        s = pd.DataFrame({"x": meta[c], "logk": logk}).dropna()
        if len(s) < 30:
            continue
        rho, pv = stats.spearmanr(s["x"], s["logk"])
        rows.append({"feature": c, "N": int(len(s)), "spearman": round(float(rho), 3), "p": float(pv)})
        print(f"  {c:18s} N={len(s):5d} Spearman={rho:+.3f} (p={pv:.2e})")
    pd.DataFrame(rows).to_csv(config.PROC_DIR / "descriptor_corr.csv", index=False)

    # ---- 图 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig_dir = config.FIG_DIR
    fig_dir.mkdir(exist_ok=True)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    d = meta.dropna(subset=["v_s_proxy"])
    ax[0].scatter(d["v_s_proxy"], np.log10(d["clarke"]), c=d["density"], cmap="viridis", s=6, alpha=0.5)
    ax[0].set_xlabel("sound-velocity proxy sqrt(B/rho) [m/s]"); ax[0].set_ylabel("log10(clarke)")
    rho, _ = stats.spearmanr(d["v_s_proxy"], np.log10(d["clarke"]))
    ax[0].set_title(f"kL vs sound velocity (Spearman={rho:+.2f}, N={len(d)})")

    d2 = meta.dropna(subset=["debye"])
    ax[1].scatter(d2["debye"], np.log10(d2["clarke"]), c=d2["bulk_vrh"], cmap="plasma", s=6, alpha=0.5)
    ax[1].set_xlabel("Debye temperature [K]"); ax[1].set_ylabel("log10(clarke)")
    rho2, _ = stats.spearmanr(d2["debye"], np.log10(d2["clarke"]))
    ax[1].set_title(f"kL vs Debye temperature (Spearman={rho2:+.2f})")

    d3 = meta.dropna(subset=["v_long"])
    ax[2].scatter(d3["v_long"], np.log10(d3["clarke"]), s=6, alpha=0.4)
    ax[2].set_xlabel("longitudinal sound velocity [m/s]"); ax[2].set_ylabel("log10(clarke)")
    rho3, _ = stats.spearmanr(d3["v_long"], np.log10(d3["clarke"]))
    ax[2].set_title(f"kL vs v_long (Spearman={rho3:+.2f})")
    plt.tight_layout(); plt.savefig(fig_dir / "kL_descriptor_scatter.png", dpi=130)
    print("saved figures/kL_descriptor_scatter.png")

    # 结构距离 vs kL 距离散点（Step 3 修复后严格对齐）
    fig2, ax = plt.subplots(figsize=(5.2, 4.4))
    d_geo = soap_distance(soap_geo); d_geo /= d_geo.max()
    d_comp = hellinger_distance(comp_frac); d_comp /= d_comp.max()
    d_struct = 0.5 * d_geo + 0.5 * d_comp
    logk_all = np.log10(meta["clarke"].values.astype(float))
    d_kL = np.abs(logk_all[:, None] - logk_all[None, :])
    n = len(meta)
    rng = np.random.RandomState(0)
    idx = rng.choice(n, size=min(4000, n), replace=False)
    sub = d_struct[np.ix_(idx, idx)]; subk = d_kL[np.ix_(idx, idx)]
    iu = np.triu_indices(len(idx), k=1)
    ax.scatter(sub[iu], subk[iu], s=1, alpha=0.15)
    rho, _ = stats.spearmanr(sub[iu], subk[iu])
    ax.set_xlabel("structure distance"); ax.set_ylabel("|log10 kL_i - log10 kL_j|")
    ax.set_title(f"structure vs kL pairwise distance (Spearman={rho:+.3f})")
    plt.tight_layout(); plt.savefig(fig_dir / "struct_vs_kL_dist.png", dpi=130)
    print("saved figures/struct_vs_kL_dist.png")


if __name__ == "__main__":
    main()
