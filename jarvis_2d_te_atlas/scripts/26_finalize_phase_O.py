# Phase O 收尾: graph QA 汇总 + 高 disagreement + spectral diagnostics + neighbor examples + 图
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import kNN_affinity

root = Path(__file__).resolve().parents[1]
K = 15

def W_of(name):
    d = np.load(root / 'graphs' / ('G_' + name + '.npz'), allow_pickle=True)
    W = d['W']
    return W.item() if hasattr(W, 'item') else W

def jids_of(name):
    return pd.read_csv(root / 'graphs' / ('G_' + name + '_nodes.csv'))['jid'].tolist()

# 1) graph QA 汇总
names = ['structure_v1','Eg_v1','electronic_n_v1','electronic_p_v1','n_transport_v1','p_transport_v1',
         'electronic_joint_sensitivity','n_transport_kappa_sensitivity','p_transport_kappa_sensitivity']
qa_rows = []
for n in names:
    W = W_of(n)
    ncomp, lab = sparse.csgraph.connected_components(W, directed=False)
    sz = np.bincount(lab)
    deg = np.asarray((W > 0).sum(axis=1)).ravel()
    qa_rows.append({'graph': 'G_' + n, 'N_nodes': W.shape[0], 'N_edges': int((W>0).nnz)//2,
                    'n_components': ncomp, 'giant_fraction': round(float(sz.max()/W.shape[0]), 4),
                    'isolated': int((deg==0).sum()), 'mean_degree': round(float(deg.mean()),2), 'k': K})
qdf = pd.DataFrame(qa_rows)
qdf.to_csv(root / 'data/audit/single_view_graph_QA.csv', index=False)
print('=== single_view_graph_QA.csv ===')
print(qdf.to_string(index=False))

# 2) 高 structure-transport disagreement (O18-O19)
d_struct = np.load(root / 'data/processed/d_struct_baseline.npy')
sj = pd.read_parquet(root / 'features/structure/geometry_soap_v1.parquet').sort_values('jid')['jid'].tolist()
smap = {j:i for i,j in enumerate(sj)}
def knn_of(D, jids, k):
    order = np.argsort(D, axis=1)[:, 1:k+1]
    return {jids[i]: set(order[i]) for i in range(len(jids))}
knn_struct = knn_of(d_struct, sj, K)
for carrier, gname in [('n','n_transport_v1'), ('p','p_transport_v1')]:
    W = W_of(gname); jids_t = jids_of(gname)
    # 从 W 重建 kNN（按 affinity 行排序 top-k）
    def knn_from_W(W, k):
        out = {}
        Wc = W.tocsr()
        for i in range(W.shape[0]):
            row = Wc[i].toarray().ravel()
            out[jids_t[i]] = set(np.argsort(-row)[:k+1].tolist()) - {i}
        return out
    knn_t = knn_from_W(W, K)
    common = sorted(set(sj) & set(jids_t))
    rows = []
    for jid in common:
        si = smap[jid]
        A = knn_struct[jid]; B = knn_t[jid]
        overlap = len(A & B) / K
        disagreement = 1 - overlap
        rows.append({'jid': jid, 'disagreement': round(disagreement, 4)})
    dd = pd.DataFrame(rows).sort_values('disagreement', ascending=False).head(50)
    dd.to_csv(root / 'data/audit' / ('high_structure_transport_disagreement_' + carrier + '.csv'), index=False)
    print(f'high_structure_transport_disagreement_{carrier}.csv: top 50 saved, max disagreement={dd.disagreement.max():.3f}')

# 3) spectral diagnostics (O2): normalized Laplacian 最小 10 个本征值
import scipy.sparse.linalg as sla
print('=== spectral diagnostics (smallest eigvals of L_sym) ===')
for n in ['structure_v1','n_transport_v1','p_transport_v1']:
    W = W_of(n)
    deg = np.asarray(W.sum(axis=1)).ravel()
    d_inv_sqrt = 1.0 / np.sqrt(np.maximum(deg, 1e-12))
    D = sparse.diags(d_inv_sqrt)
    L = sparse.identity(W.shape[0]) - D @ W @ D
    try:
        ev = sla.eigsh(L, k=8, which='SM', return_eigenvectors=False)
        ev = np.sort(ev)
        print(f'  {n}: smallest 8 eigvals = {np.round(ev, 4).tolist()}')
    except Exception as e:
        print(f'  {n}: eigsh failed ({e})')

# 4) neighbor examples (O20): 每个主 view 抽 3 个 anchor, 列 10 近邻
formula = pd.read_parquet(root / 'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
rng = np.random.RandomState(0)
for gname, outname in [('structure_v1','structure'), ('electronic_n_v1','electronic_n'), ('electronic_p_v1','electronic_p'), ('n_transport_v1','transport_n'), ('p_transport_v1','transport_p')]:
    W = W_of(gname); jids = jids_of(gname)
    lines = []
    anchors = rng.choice(len(jids), 3, replace=False)
    for ai in anchors:
        row = W[ai].toarray().ravel()
        top = np.argsort(-row)[:10]
        lines.append('## anchor ' + jids[ai] + ' (' + str(formula.get(jids[ai],'')) + ')')
        for t in top:
            lines.append('  - ' + jids[t] + ' (' + str(formula.get(jids[t],'')) + ') w=%.3f' % row[t])
        lines.append('')
    (root / 'reports' / ('neighbor_examples_' + outname + '.md')).write_text('\n'.join(lines), encoding='utf-8')
print('wrote neighbor_examples_*.md (5 files)')
