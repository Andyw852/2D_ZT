"""kappa_L verify: cross-view (structure/electronic vs kL) neighbor-overlap + distance correlation.
Mirrors 24_cross_view_analysis.py, structure view matches scripts 20/22 (geo-SOAP + composition Hellinger).
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import RobustScaler
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import hellinger_distance, soap_distance

root = Path(__file__).resolve().parents[1]
K = 10
rng = np.random.RandomState(42)
N_PERM = 2000

df = pd.read_parquet(root / "features" / "kl_verify" / "kl_views.parquet")
soap_mean = np.load(root / "data" / "processed" / "kl_soap_geo.npy")
frac = np.load(root / "data" / "processed" / "kl_comp_frac.npy")
assert len(df) == soap_mean.shape[0] == frac.shape[0]
print("N materials:", len(df))

# ---- distances ----
d_geo = soap_distance(soap_mean)
d_comp = hellinger_distance(frac)
d_geo_n = d_geo / d_geo.max()
d_comp_n = d_comp / d_comp.max()
d_struct = 0.5 * d_geo_n + 0.5 * d_comp_n

def D_from_cols(cols, log=False, scale=True):
    sub = df.dropna(subset=cols).reset_index(drop=True)
    X = sub[cols].values.astype(float)
    if log:
        X = np.log10(X)
    if scale and X.shape[1] > 1:
        X = RobustScaler().fit_transform(X)
    return squareform(pdist(X)), sub["jid"].tolist()

views = {
    "Structure": (d_struct, df["jid"].tolist()),
    "Structure_geo": (d_geo_n, df["jid"].tolist()),
    "Structure_comp": (d_comp_n, df["jid"].tolist()),
    "Eg": D_from_cols(["Eg_opt"], scale=False),
    "Electronic": D_from_cols(["Eg_opt", "m_elec", "m_hole"]),
    "Elec_Eg+me": D_from_cols(["Eg_opt", "m_elec"]),
    "Elec_Eg+mh": D_from_cols(["Eg_opt", "m_hole"]),
    "kL": D_from_cols(["kL_300"], log=True, scale=False),
    "Elastic(B,G)": D_from_cols(["B_kv", "G_gv"]),
}

pairs = [
    ("Structure", "kL"),
    ("Structure_geo", "kL"),
    ("Structure_comp", "kL"),
    ("Eg", "kL"),
    ("Electronic", "kL"),
    ("Elec_Eg+me", "kL"),
    ("Elec_Eg+mh", "kL"),
    ("Elastic(B,G)", "kL"),   # physical positive control (sound velocity / Debye proxy)
    ("Structure", "Eg"),       # sanity
    ("Structure", "Electronic"),
    ("Electronic", "Elastic(B,G)"),
]

def knn_sets_by_jid(D, jids, k):
    """Return neighbor identities, not local row numbers.

    Local row indices are incomparable when two views have different coverage.
    The previous implementation intersected those indices and therefore biased
    every overlap involving the 82-material electronic/elastic subsets.
    """
    out = {}
    for i, jid in enumerate(jids):
        order = np.argsort(D[i])
        order = order[order != i][:k]
        out[jid] = {jids[z] for z in order}
    return out

knn = {name: knn_sets_by_jid(D, jids, K) for name, (D, jids) in views.items()}
jmap = {name: {j: i for i, j in enumerate(jids)} for name, (_, jids) in views.items()}

print(f"\n=== cross-view neighbor overlap (k={K}) + permutation baseline ({N_PERM}x) ===")
rows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    if len(common) < 20:
        print(f"  {A} vs {B}: N={len(common)} <20 skip")
        continue
    ov = np.mean([len(knn[A][j] & knn[B][j]) / K for j in common])
    null = []
    for _ in range(N_PERM):
        perm = rng.permutation(len(common))
        ovn = np.mean([len(knn[A][common[i]] & knn[B][common[perm[i]]]) / K for i in range(len(common))])
        null.append(ovn)
    null = np.array(null)
    z = (ov - null.mean()) / (null.std() + 1e-12)
    p = (1 + (null >= ov).sum()) / (N_PERM + 1)
    rows.append({"pair": f"{A} vs {B}", "N": len(common), "overlap": round(float(ov),4),
                 "null_mean": round(float(null.mean()),4), "z": round(float(z),2), "p": round(float(p),4)})
    print(f"  {A:18s} vs {B:6s}: N={len(common):3d} overlap={ov:.4f} null={null.mean():.4f}+-{null.std():.4f} z={z:5.1f} p={p:.4f}")

pd.DataFrame(rows).to_csv(root / "data" / "audit" / "kl_view_overlap.csv", index=False)

print("\n=== cross-view distance correlation (Spearman + node-permutation Mantel test) ===")
drows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    if len(common) < 20: continue
    ia = [jmap[A][j] for j in common]; ib = [jmap[B][j] for j in common]
    DA = views[A][0][np.ix_(ia, ia)]; DB = views[B][0][np.ix_(ib, ib)]
    n = len(common); iu = np.triu_indices(n, k=1)
    rho, naive_p = stats.spearmanr(DA[iu], DB[iu])
    null = []
    for _ in range(N_PERM):
        perm = rng.permutation(n)
        DBp = DB[np.ix_(perm, perm)]
        null.append(stats.spearmanr(DA[iu], DBp[iu]).statistic)
    null = np.asarray(null)
    mantel_p = (1 + (np.abs(null) >= abs(rho)).sum()) / (N_PERM + 1)
    drows.append({"pair": f"{A} vs {B}", "N": n, "spearman": round(float(rho),4),
                  "mantel_p": float(mantel_p), "naive_pairwise_p_invalid": float(naive_p)})
    print(f"  {A:18s} vs {B:6s}: N={n:3d} Spearman={rho:+.4f} Mantel-p={mantel_p:.4g}")

pd.DataFrame(drows).to_csv(root / "data" / "audit" / "kl_view_distance_corr.csv", index=False)

print("\n=== kL distribution (W/mK) ===")
print(df["kL_300"].describe().round(3).to_string())
print("log10(kL) range:", round(float(np.log10(df["kL_300"].min())),2), "-", round(float(np.log10(df["kL_300"].max())),2))
