# Phase O: PF/Eg external label diagnostics + leave-one-out + structure ablation + near-duplicates
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import squareform, pdist
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import hellinger_distance, soap_distance

root = Path(__file__).resolve().parents[1]
rng = np.random.RandomState(42)
K = 15

def smoothness(W, y):
    y = np.asarray(y, dtype=float)
    ii, jj = W.nonzero()
    w = np.asarray(W[ii, jj]).ravel()
    dy2 = (y[ii] - y[jj]) ** 2
    ok = np.isfinite(dy2)
    return float((w[ok] * dy2[ok]).sum() / w[ok].sum()) if ok.any() else np.nan

def null_smoothness(W, y, n=1000):
    obs = smoothness(W, y)
    null = np.array([smoothness(W, rng.permutation(y)) for _ in range(n)])
    z = (obs - null.mean()) / (null.std() + 1e-12)
    p = (null <= obs).mean()
    return obs, null.mean(), null.std(), z, p

def W_from_npz(name):
    d = np.load(root / 'graphs' / ('G_' + name + '.npz'), allow_pickle=True)
    W = d['W']
    if hasattr(W, 'item'):
        W = W.item()
    return W

def jids_of(name):
    return pd.read_csv(root / 'graphs' / ('G_' + name + '_nodes.csv'))['jid'].tolist()

pf_n = pd.read_parquet(root / 'features/transport/n_transport_features_v1.parquet')[['jid','PF_mean']].rename(columns={'PF_mean':'PF'})
pf_p = pd.read_parquet(root / 'features/transport/p_transport_features_v1.parquet')[['jid','PF_mean']].rename(columns={'PF_mean':'PF'})
edf = pd.read_parquet(root / 'features/electronic/electronic_features_v1.parquet')

print('=== PF smoothness (log10(PF+eps)) + null test ===')
eps = 1e-3
for pf, gname in [(pf_n, 'n_transport_v1'), (pf_p, 'p_transport_v1')]:
    W = W_from_npz(gname); jids = jids_of(gname)
    sub = pf.set_index('jid').reindex(jids)['PF'].values
    y = np.log10(np.maximum(sub, eps))
    obs, nm, ns, z, p = null_smoothness(W, y)
    print('  %s: smoothness=%.4f null=%.4f+-%.4f z=%.2f p=%.4f' % (gname, obs, nm, ns, z, p))

print('=== PF smoothness across views (n-type) ===')
for gname, view in [('structure_v1','Structure'), ('electronic_n_v1','Electronic_n'), ('n_transport_v1','Transport_n')]:
    W = W_from_npz(gname); jids = jids_of(gname)
    sub = pf_n.set_index('jid').reindex(jids)['PF'].values
    y = np.log10(np.maximum(sub, eps))
    obs, nm, ns, z, p = null_smoothness(W, y)
    print('  %s: smoothness=%.4f z=%.2f p=%.4f' % (view, obs, z, p))

print('=== Eg smoothness on structure graph ===')
W = W_from_npz('structure_v1'); jids = jids_of('structure_v1')
egv = edf.set_index('jid').reindex(jids)['Eg_optb88vdw'].values
obs, nm, ns, z, p = null_smoothness(W, egv)
print('  Structure graph Eg smoothness=%.4f z=%.2f p=%.4f' % (obs, z, p))

print('=== n-transport feature leave-one-out (k=15) ===')
V1 = ['S_median','S_MAD','S_sign_fraction','log_sigma_dom_geo','D_sigma','A_sigma_dom']
tdf = pd.read_parquet(root / 'features/transport/n_transport_features_v1.parquet')
def dist_of(df, cols):
    sub = df.dropna(subset=cols).reset_index(drop=True)
    X = RobustScaler().fit_transform(sub[cols].values)
    return squareform(pdist(X)), sub['jid'].tolist()
Dfull, jids = dist_of(tdf, V1)
def knn_overlap(da, db, k):
    ka = np.argsort(da, axis=1)[:, 1:k+1]; kb = np.argsort(db, axis=1)[:, 1:k+1]
    return np.mean([len(set(ka[i]) & set(kb[i])) / k for i in range(len(da))])
loo = []
for f in V1:
    cols = [c for c in V1 if c != f]
    Dl, _ = dist_of(tdf, cols)
    ov = knn_overlap(Dfull, Dl, K)
    loo.append({'feature': f, 'knn_overlap_when_removed': round(ov, 4)})
    print('  remove %s: kNN overlap = %.4f' % (f, ov))
pd.DataFrame(loo).to_csv(root / 'data/audit/transport_feature_leave_one_out.csv', index=False)

print('=== structure block ablation ===')
gdf = pd.read_parquet(root / 'features/structure/geometry_soap_v1.parquet').sort_values('jid')
d_geo = soap_distance(gdf.filter(regex='soap6_mean_').values)
F = np.array([json.loads(x) for x in pd.read_parquet(root / 'features/structure/composition_fraction.parquet').sort_values('jid')['fraction']])
d_comp = hellinger_distance(F)
d_geo_n = d_geo / d_geo.max(); d_comp_n = d_comp / d_comp.max()
d_comb = 0.5 * d_geo_n + 0.5 * d_comp_n
print('  geo vs comp kNN(20) = %.4f' % knn_overlap(d_geo, d_comp, 20))
print('  geo vs combined kNN(20) = %.4f' % knn_overlap(d_geo, d_comb, 20))
print('  comp vs combined kNN(20) = %.4f' % knn_overlap(d_comp, d_comb, 20))

print('=== structure near-duplicates ===')
soap_df = gdf.reset_index(drop=True); jids_all = soap_df['jid'].tolist()
formula = pd.read_parquet(root / 'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
iu = np.triu_indices(len(jids_all), k=1)
dists = d_comb[iu]
idx = np.argsort(dists)[:100]
nd_rows = []
for t in idx:
    i, j = iu[0][t], iu[1][t]
    nd_rows.append({'jid_A': jids_all[i], 'jid_B': jids_all[j], 'formula_A': formula[jids_all[i]], 'formula_B': formula[jids_all[j]], 'd_struct': round(float(dists[t]), 6), 'same_formula': formula[jids_all[i]] == formula[jids_all[j]]})
ndf = pd.DataFrame(nd_rows)
ndf.to_csv(root / 'data/audit/structure_near_duplicates.csv', index=False)
print(ndf.head(20).to_string(index=False))