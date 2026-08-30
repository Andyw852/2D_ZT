"""Phase O: cross-view neighbor overlap + distance correlation + random baseline + PF/Eg diagnostics。"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import cdist, pdist, squareform
sys.path.insert(0, str(Path(__file__).resolve().parent))

root = Path(__file__).resolve().parents[1]
K = 15
rng = np.random.RandomState(42)

def scaled_dist_from(df, cols):
    sub = df.dropna(subset=cols).reset_index(drop=True)
    X = RobustScaler().fit_transform(sub[cols].values)
    return squareform(pdist(X)), sub["jid"].tolist()

# ---- 载入各 view 距离 + jids ----
views = {}
# Structure
d_struct = np.load(root / "data" / "processed" / "d_struct_baseline.npy")
soap_df = pd.read_parquet(root / "features" / "structure" / "geometry_soap_v1.parquet").sort_values("jid")
views["Structure"] = (d_struct, soap_df["jid"].tolist())
# Eg
edf = pd.read_parquet(root / "features" / "electronic" / "electronic_features_v1.parquet")
eg = edf[["jid","Eg_optb88vdw"]].dropna()
views["Eg"] = (squareform(pdist(eg[["Eg_optb88vdw"]].values)), eg["jid"].tolist())
# Electronic-n/p
for name, cols in [("Electronic_n", ["Eg_optb88vdw","m_elec_median"]),
                   ("Electronic_p", ["Eg_optb88vdw","m_hole_median"])]:
    D, jids = scaled_dist_from(edf, cols)
    views[name] = (D, jids)
# Transport n/p
V1 = ["S_median","S_MAD","S_sign_fraction","log_sigma_dom_geo","D_sigma","A_sigma_dom"]
for name, f in [("Transport_n", "n_transport_features_v1.parquet"),
                ("Transport_p", "p_transport_features_v1.parquet")]:
    tdf = pd.read_parquet(root / "features" / "transport" / f)
    D, jids = scaled_dist_from(tdf, V1)
    views[name] = (D, jids)

def knn_sets(D, k):
    n = D.shape[0]
    order = np.argsort(D, axis=1)[:, 1:k+1]
    return [set(order[i]) for i in range(n)]

knn = {name: knn_sets(D, K) for name, (D, _) in views.items()}
jmap = {name: {j: i for i, j in enumerate(jids)} for name, (_, jids) in views.items()}

pairs = [("Structure","Eg"),("Structure","Electronic_n"),("Structure","Electronic_p"),
         ("Structure","Transport_n"),("Structure","Transport_p"),
         ("Electronic_n","Transport_n"),("Electronic_p","Transport_p"),
         ("Transport_n","Transport_p")]

print("=== cross-view neighbor overlap (k=15) + random baseline ===")
rows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    ia = [jmap[A][j] for j in common]; ib = [jmap[B][j] for j in common]
    ov = np.mean([len(knn[A][ia[i]] & knn[B][ib[i]]) / K for i in range(len(common))])
    # random baseline: 打乱 B 的 JID 对应 1000 次
    null = []
    order_b = np.array(ib)
    for _ in range(1000):
        perm = rng.permutation(len(common))
        ovn = np.mean([len(knn[A][ia[i]] & knn[B][order_b[perm[i]]]) / K for i in range(len(common))])
        null.append(ovn)
    null = np.array(null)
    z = (ov - null.mean()) / (null.std() + 1e-12)
    p = (null >= ov).mean()
    rows.append({"pair": f"{A} vs {B}", "N_common": len(common), "overlap": round(ov,4),
                 "null_mean": round(null.mean(),4), "null_std": round(null.std(),4),
                 "z_score": round(z,2), "empirical_p": round(p,4)})
    print(f"  {A} vs {B}: N={len(common)} overlap={ov:.4f} null={null.mean():.4f}+-{null.std():.4f} z={z:.1f} p={p:.4f}")
ovdf = pd.DataFrame(rows)
ovdf.to_csv(root / "data" / "audit" / "view_neighbor_overlap.csv", index=False)

# ---- distance correlation ----
print("\n=== cross-view distance correlation (Spearman, sampled) ===")
drows = []
for A, B in pairs:
    common = sorted(set(views[A][1]) & set(views[B][1]))
    ia = [jmap[A][j] for j in common]; ib = [jmap[B][j] for j in common]
    DA = views[A][0][np.ix_(ia, ia)]; DB = views[B][0][np.ix_(ib, ib)]
    n = len(common); iu = np.triu_indices(n, k=1)
    da = DA[iu]; db = DB[iu]
    if len(da) > 200000:
        idx = rng.choice(len(da), 200000, replace=False)
        da = da[idx]; db = db[idx]
    rho = stats.spearmanr(da, db)[0]
    drows.append({"pair": f"{A} vs {B}", "N_common": n, "spearman": round(rho,4)})
    print(f"  {A} vs {B}: N={n} Spearman={rho:.4f}")
dd = pd.DataFrame(drows)
dd.to_csv(root / "data" / "audit" / "view_distance_correlation.csv", index=False)

# ---- view similarity matrices (heatmap) ----
view_names = ["Structure","Eg","Electronic_n","Electronic_p","Transport_n","Transport_p"]
OM = np.zeros((6,6)); DM = np.zeros((6,6))
for i, A in enumerate(view_names):
    for j, B in enumerate(view_names):
        common = sorted(set(views[A][1]) & set(views[B][1]))
        ia = [jmap[A][j] for j in common]; ib = [jmap[B][j] for j in common]
        ov = np.mean([len(knn[A][ia[t]] & knn[B][ib[t]]) / K for t in range(len(common))]) if common else np.nan
        OM[i,j] = ov
        DA = views[A][0][np.ix_(ia, ia)]; DB = views[B][0][np.ix_(ib, ib)]
        n = len(common); iu = np.triu_indices(n, k=1)
        if n > 1:
            DM[i,j] = stats.spearmanr(DA[iu], DB[iu])[0]
        else:
            DM[i,j] = np.nan
np.save(root / "data" / "processed" / "view_overlap_matrix.npy", OM)
np.save(root / "data" / "processed" / "view_distance_corr_matrix.npy", DM)

import pandas as _pd
print("view neighbor-overlap matrix:")
print(_pd.DataFrame(OM, index=view_names, columns=view_names).round(3).to_string())
print("view distance-correlation matrix:")
print(_pd.DataFrame(DM, index=view_names, columns=view_names).round(3).to_string())
