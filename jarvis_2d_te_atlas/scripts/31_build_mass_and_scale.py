"""Phase P-A/P-B: mass-only graphs + layer strength scaling。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from datetime import datetime, timezone
sys.path.insert(0, str(Path(__file__).resolve().parent))
from multiview_utils import kNN_affinity_from_dist, scale_layer, load_W, jids_of
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import squareform, pdist
from scipy.sparse.csgraph import connected_components

root = Path(__file__).resolve().parents[1]
gdir = root / 'graphs'
now = datetime.now(timezone.utc).isoformat()
edf = pd.read_parquet(root / 'features/electronic/electronic_features_v1.parquet')

# P-A: mass-only graphs
for name, col in [('m_electron_v1','m_elec_median'), ('m_hole_v1','m_hole_median')]:
    sub = edf.dropna(subset=[col]).reset_index(drop=True)
    X = RobustScaler().fit_transform(sub[[col]].values)
    D = squareform(pdist(X))
    W = kNN_affinity_from_dist(D, 15)
    np.savez_compressed(gdir / ('G_' + name + '.npz'), W=W)
    pd.DataFrame({'jid': sub['jid'], 'degree': np.asarray((W>0).sum(axis=1)).ravel()}).to_csv(gdir / ('G_' + name + '_nodes.csv'), index=False)
    ncomp, lab = connected_components(W, directed=False)
    sz = np.bincount(lab)
    meta = {'graph_name':'G_'+name,'view':'mass','features':[col],'feature_transform':'RobustScaler','distance_metric':'Euclidean','k':15,
            'kernel':'local-scale Gaussian','symmetrization':'union','normalization':'none','N_nodes':W.shape[0],
            'n_components':ncomp,'giant_component_fraction':round(float(sz.max()/W.shape[0]),4),'creation_date':now,'random_seed':42}
    (gdir / ('G_' + name + '_metadata.json')).write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print(f'G_{name}: N={W.shape[0]} giant={sz.max()/W.shape[0]:.4f} comp={ncomp}')

# P-B: layer scaling
print('\n=== layer scaling (mean strength -> 1) ===')
layers_scale = [('structure_v1',1.0), ('Eg_v1',1.0), ('m_electron_v1',1.0), ('m_hole_v1',1.0),
                ('n_transport_v1',1.0), ('p_transport_v1',1.0)]
rows = []
scaled = {}
for name, alpha in layers_scale:
    W = load_W(name, gdir)
    Ws, ms = scale_layer(W)
    scaled[name] = Ws
    deg = np.asarray(W.sum(axis=1)).ravel()
    degs = np.asarray(Ws.sum(axis=1)).ravel()
    rows.append({'view': name, 'N_nodes': W.shape[0], 'N_edges': int((W>0).nnz)//2,
                 'mean_strength_before': round(float(deg.mean()),4), 'median_strength_before': round(float(np.median(deg)),4),
                 'mean_strength_after': round(float(degs.mean()),4), 'scale_factor': round(1.0/ms,6)})
sdf = pd.DataFrame(rows)
sdf.to_csv(root / 'data/audit/multiview_layer_scale.csv', index=False)
print(sdf.to_string(index=False))
# 保存 scaled graphs
np.savez_compressed(gdir / 'multiview_scaled.npz', **{k: v for k, v in scaled.items()})
print('saved multiview_scaled.npz')
