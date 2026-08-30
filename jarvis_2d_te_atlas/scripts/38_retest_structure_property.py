"""重测 Structure→property 相关性（更强描述符 Magpie+SOAP vs 旧 Hellinger+SOAP）。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import squareform, pdist
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import hellinger_distance, soap_distance

root = Path(__file__).resolve().parents[1]
gdf = pd.read_parquet(root/'features/structure/geometry_soap_v1.parquet').sort_values('jid').reset_index(drop=True)
mdf = pd.read_parquet(root/'features/structure/composition_magpie.parquet').sort_values('jid').reset_index(drop=True)
cdf = pd.read_parquet(root/'features/structure/composition_fraction.parquet').sort_values('jid').reset_index(drop=True)
jids = gdf['jid'].tolist()

# 距离
d_geo = soap_distance(gdf.filter(regex='soap6_mean_').values)
d_comp_old = hellinger_distance(np.array([json.loads(x) for x in cdf['fraction']]))
Xm = mdf.iloc[:, 1:].values
Xm = np.nan_to_num(Xm)
Xm = RobustScaler().fit_transform(Xm)
d_comp_new = squareform(pdist(Xm))  # Magpie 成分距离

def norm(D): return D / (D.max() + 1e-12)
d_geo_n = norm(d_geo); d_comp_old_n = norm(d_comp_old); d_comp_new_n = norm(d_comp_new)

# 属性
tdf_n = pd.read_parquet(root/'features/transport/n_transport_features_v1.parquet').set_index('jid')
edf = pd.read_parquet(root/'features/electronic/electronic_features_v1.parquet').set_index('jid')
props = {
    'Eg': edf['Eg_optb88vdw'],
    'm_elec': edf['m_elec_median'],
    'S_median': tdf_n['S_median'],
    'log_sigma_dom': tdf_n['log_sigma_dom_geo'],
    'PF_mean': tdf_n['PF_mean'],
}

def dist_corr(D, y, jids):
    # 在共同非缺失 jid 上，比较 pairwise D 与 pairwise |dy|
    sub = y.dropna()
    common = [j for j in jids if j in sub.index]
    idx = [jids.index(j) for j in common]
    Dsub = D[np.ix_(idx, idx)]
    yv = sub.reindex(common).values
    iu = np.triu_indices(len(common), k=1)
    d = Dsub[iu]
    dy = np.abs(yv[:,None] - yv[None,:])[iu]
    r = stats.spearmanr(d, dy)[0]
    return r, len(common)

print('=== Structure→property distance correlation (Spearman) ===')
print(f'{"property":<16} {"old(Hellinger+SOAP)":<22} {"new(Magpie+SOAP)":<22} {"N"}')
for pname, y in props.items():
    r_old, n = dist_corr(0.5*d_geo_n + 0.5*d_comp_old_n, y, jids)
    r_new, n = dist_corr(0.5*d_geo_n + 0.5*d_comp_new_n, y, jids)
    print(f'{pname:<16} {r_old:<22.4f} {r_new:<22.4f} {n}')

# 旧 vs 新 结构距离的一致性
iu = np.triu_indices(len(jids), k=1)
d_old_struct = 0.5*d_geo_n + 0.5*d_comp_old_n
d_new_struct = 0.5*d_geo_n + 0.5*d_comp_new_n
rho = stats.spearmanr(d_old_struct[iu], d_new_struct[iu])[0]
print(f'\nold vs new structure distance Spearman = {rho:.4f}')

# 单独看 composition 距离（旧 Hellinger vs 新 Magpie）对 property 的相关
print('\n=== 仅 composition 距离对 property 的相关 ===')
for pname, y in props.items():
    r_old, _ = dist_corr(d_comp_old_n, y, jids)
    r_new, _ = dist_corr(d_comp_new_n, y, jids)
    print(f'{pname:<16} Hellinger={r_old:.4f}  Magpie={r_new:.4f}')
