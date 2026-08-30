"""MP 数据: 跨视图(结构/弹性/电子 vs 晶格热导率) 近邻重叠 + 距离相关。
方法论与 jarvis_2d_te_atlas/scripts/kl_verify_03_crossview.py 一致。
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import RobustScaler

root = Path(__file__).resolve().parents[1]
mp = root / "mp_kappaL"
sys.path.insert(0, str(root / "jarvis_2d_te_atlas" / "scripts"))
from graph_utils import hellinger_distance, soap_distance

K = 10
rng = np.random.RandomState(42)
N_PERM = 300
N_PERM_MANTEL = 100

meta = pd.read_parquet(mp / "processed" / "views_meta.parquet")
soap_geo = np.load(mp / "processed" / "soap_geo.npy").astype(np.float64)
comp_frac = np.load(mp / "processed" / "comp_frac.npy").astype(np.float64)
elec = pd.read_parquet(mp / "processed" / "electronic_jarvis.parquet")
meta = meta.merge(elec, on="material_id", how="left")

print("N materials (structure/kL):", len(meta))

# ---- 距离矩阵 ----
d_geo = soap_distance(soap_geo).astype(np.float32)
d_comp = hellinger_distance(comp_frac).astype(np.float32)
d_geo_n = d_geo / d_geo.max()
d_comp_n = d_comp / d_comp.max()
d_struct = (0.5 * d_geo_n + 0.5 * d_comp_n).astype(np.float32)

ids = meta["material_id"].tolist()

def D_from_cols(cols, log=False, scale=True, filt=None):
    sub = meta.copy()
    if filt is not None:
        sub = sub[filt(sub)]
    sub = sub.dropna(subset=cols).reset_index(drop=True)
    X = sub[cols].values.astype(float)
    if log:
        X = np.log10(X)
    if scale and X.shape[1] > 1:
        X = RobustScaler().fit_transform(X)
    return squareform(pdist(X)).astype(np.float32), sub["material_id"].tolist()

# 弹性视图: 过滤异常体/剪切模量 (物理合理范围)
elastic_filt = lambda s: (s["bulk_vrh"] > 0) & (s["bulk_vrh"] < 1000) & (s["shear_vrh"] > 0) & (s["shear_vrh"] < 1000)
# snyder_acoustic 过滤数值爆炸 (>1e4 W/mK 视为拟合失败)
snyder_filt = lambda s: (s["snyder_acoustic"] > 0) & (s["snyder_acoustic"] < 1e4)

views = {
    "Structure":      (d_struct, ids),
    "Structure_geo":  (d_geo_n, ids),
    "Structure_comp": (d_comp_n, ids),
    "kL_clarke":      D_from_cols(["clarke"], log=True, scale=False),
    "kL_cahill":      D_from_cols(["cahill"], log=True, scale=False),
    "kL_snyder":      D_from_cols(["snyder_acoustic"], log=True, scale=False, filt=snyder_filt),
    "Elastic":        D_from_cols(["bulk_vrh", "shear_vrh", "debye"], filt=elastic_filt),
    "Eg":             D_from_cols(["gap_opt"], scale=False),
    "Electronic":     D_from_cols(["gap_opt", "m_elec", "m_hole"]),
    "Eg+me":          D_from_cols(["gap_opt", "m_elec"]),
    "Eg+mh":          D_from_cols(["gap_opt", "m_hole"]),
}

pairs = [
    ("Structure", "kL_clarke"),
    ("Structure_geo", "kL_clarke"),
    ("Structure_comp", "kL_clarke"),
    ("Elastic", "kL_clarke"),
    ("Eg", "kL_clarke"),
    ("Electronic", "kL_clarke"),
    ("Structure", "kL_snyder"),
    ("Elastic", "kL_snyder"),
    ("Eg", "kL_snyder"),
    ("Structure", "Eg"),
    ("Structure", "Electronic"),
    ("Electronic", "Elastic"),
    ("Elastic", "Structure"),
    ("kL_clarke", "kL_cahill"),
    ("kL_clarke", "kL_snyder"),
]

def knn_sets(D, jids, k):
    out = {}
    n = D.shape[0]
    for i in range(n):
        row = D[i]
        idx = np.argpartition(row, k + 1)[: k + 1]
        nn = idx[idx != i][:k]
        out[jids[i]] = {jids[z] for z in nn}
    return out

print("computing kNN sets ...")
knn = {name: knn_sets(D, jids, K) for name, (D, jids) in views.items()}
jmap = {name: {j: i for i, j in enumerate(jids)} for name, (_, jids) in views.items()}

print(f"=== cross-view neighbor overlap (k={K}) + permutation ({N_PERM}x) ===")
rows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    if len(common) < 20:
        print(f"  {A} vs {B}: N={len(common)} <20 skip"); continue
    ov = np.mean([len(knn[A][j] & knn[B][j]) / K for j in common])
    null = []
    for _ in range(N_PERM):
        perm = rng.permutation(len(common))
        ovn = np.mean([len(knn[A][common[i]] & knn[B][common[perm[i]]]) / K for i in range(len(common))])
        null.append(ovn)
    null = np.array(null)
    z = (ov - null.mean()) / (null.std() + 1e-12)
    p = (1 + (null >= ov).sum()) / (N_PERM + 1)
    enrich = ov / (null.mean() + 1e-12)
    rows.append({"pair": f"{A} vs {B}", "N": len(common), "overlap": round(float(ov),4),
                 "null_mean": round(float(null.mean()),4), "z": round(float(z),1),
                 "p": round(float(p),4), "enrichment": round(float(enrich),2)})
    print(f"  {A:16s} vs {B:12s}: N={len(common):5d} overlap={ov:.4f} null={null.mean():.4f} z={z:6.1f} p={p:.4f} enrich={enrich:.2f}x")

pd.DataFrame(rows).to_csv(mp / "processed" / "view_overlap.csv", index=False)

print(f"=== cross-view distance correlation (Spearman on sampled pairs + Mantel {N_PERM_MANTEL}x) ===")
headline = {("Structure","kL_clarke"), ("Elastic","kL_clarke"), ("Eg","kL_clarke"), ("Electronic","kL_clarke")}
M_SAMP = 1_000_000
drows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    if len(common) < 20: continue
    ia = [jmap[A][j] for j in common]; ib = [jmap[B][j] for j in common]
    DA = views[A][0][np.ix_(ia, ia)]; DB = views[B][0][np.ix_(ib, ib)]
    n = len(common)
    iu = np.triu_indices(n, k=1)
    nsamp = min(M_SAMP, len(iu[0]))
    sel = rng.choice(len(iu[0]), size=nsamp, replace=False)
    ip, jp = iu[0][sel], iu[1][sel]
    da = DA[ip, jp]; db = DB[ip, jp]
    rho, _ = stats.spearmanr(da, db)
    mantel_p = np.nan
    if (A, B) in headline and n <= 13000:
        null = []
        for _ in range(N_PERM_MANTEL):
            perm = rng.permutation(n)
            DBp = DB[np.ix_(perm, perm)]
            null.append(stats.spearmanr(da, DBp[ip, jp]).statistic)
        null = np.asarray(null)
        mantel_p = (1 + (np.abs(null) >= abs(rho)).sum()) / (N_PERM_MANTEL + 1)
    drows.append({"pair": f"{A} vs {B}", "N": n, "spearman": round(float(rho),4), "mantel_p": float(mantel_p)})
    print(f"  {A:16s} vs {B:12s}: N={n:5d} Spearman={rho:+.4f} Mantel-p={mantel_p:.4g}")

pd.DataFrame(drows).to_csv(mp / "processed" / "view_distance_corr.csv", index=False)

print("=== kL distribution (clarke, W/mK) ===")
print(meta["clarke"].describe().round(3).to_string())
print("log10(clarke) range:", round(float(np.log10(meta['clarke'].min())),2), "-", round(float(np.log10(meta['clarke'].max())),2))
