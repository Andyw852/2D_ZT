"""Phase R-S: preservation real vs random + joint tension + PF smoothness + metal/semi。"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse, stats
sys.path.insert(0, 'scripts')
from multiview_utils import build_supra, joint_embedding, load_W, jids_of

root = Path('.')
gdir = root / 'graphs'
N_EIGS = 20; LAM = 0.3; K = 15
rng = np.random.RandomState(42)
d = np.load(gdir / 'multiview_scaled.npz', allow_pickle=True)
scaled = {k: d[k].item() if hasattr(d[k], 'item') else d[k] for k in d.files}
jids_map = {n: jids_of(n, gdir) for n in ['structure_v1','Eg_v1','m_electron_v1','m_hole_v1','n_transport_v1','p_transport_v1']}
all_jids = sorted(jids_map['structure_v1'])

def layers_for(cfg):
    m = {'Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1','n_transport':'n_transport_v1','p_transport':'p_transport_v1'}
    return [(m[v], scaled[m[v]], jids_map[m[v]]) for v in cfg]

def embed_and_eval(layers, lam, shuffle=False):
    layers2 = layers
    if shuffle:
        layers2 = []
        for name, W, jids in layers:
            perm = rng.permutation(len(jids)); jids2 = [jids[p] for p in perm]
            layers2.append((name, W, jids2))
    A, offsets, cid = build_supra(layers2, all_jids, lam)
    ev, V = joint_embedding(A, N_EIGS)
    if V is None: return None
    coords = V[:, 1:N_EIGS+1]
    cons = coords[:len(all_jids)]
    # transport preservation
    tname = 'n_transport_v1' if any('n_transport' in x[0] for x in layers) else 'p_transport_v1'
    for name, W, jids in layers2:
        if name == tname:
            off = offsets[name]
            vc = coords[off:off+len(jids)]
            Dv = np.linalg.norm(vc[:,None,:]-vc[None,:,:], axis=-1)
            order = np.argsort(Dv, axis=1)[:, 1:K+1]
            joint_nn = {jids[i]: set(order[i].tolist()) for i in range(len(jids))}
            Wc = W.tocsr()
            orig_nn = {}
            for i in range(W.shape[0]):
                row = Wc[i].toarray().ravel()
                orig_nn[jids[i]] = set(np.argsort(-row)[:K+1].tolist()) - {i}
            P = np.mean([len(orig_nn[j] & joint_nn[j])/K for j in jids])
    # tension (absolute)
    ts = []
    for name, W, jids in layers2:
        jmap = {j:i for i,j in enumerate(jids)}; off = offsets[name]
        for jid in jids:
            if jid in cid: ts.append(np.linalg.norm(cons[cid[jid]] - coords[off+jmap[jid]]))
    return {'P_transport': P, 'median_tension': float(np.median(ts)), 'coords': coords, 'offsets': offsets, 'cid': cid, 'layers': layers2}

print('=== preservation real vs random (transport) ===', flush=True)
rows = []
for carrier, cfg in [('n',['Eg','m_electron','n_transport']), ('p',['Eg','m_hole','p_transport'])]:
    layers = layers_for(cfg)
    real = embed_and_eval(layers, LAM, False)
    null = [embed_and_eval(layers, LAM, True)['P_transport'] for _ in range(30)]
    null = np.array(null)
    z = (real['P_transport'] - null.mean())/(null.std()+1e-12)
    p = (null >= real['P_transport']).mean()
    rows.append({'carrier': carrier, 'P_real': round(real['P_transport'],4), 'P_null_mean': round(float(null.mean()),4),
                 'P_null_std': round(float(null.std()),4), 'z': round(float(z),2), 'p': round(float(p),4)})
    print(f"  {carrier}: P_transport real={real['P_transport']:.4f} vs null={null.mean():.4f}+-{null.std():.4f} z={z:.1f} p={p:.4f}", flush=True)
pd.DataFrame(rows).to_csv(root/'data/audit/random_anchor_preservation.csv', index=False)

# PF smoothness on joint（一致性检查，非独立验证）：
# PF_mean = S²σ 由 nseeb×ncond 定义，而 transport 图特征本身嵌入了 S 与 σ，
# 故 PF 在联合流形上平滑是「同源一致性」，不构成对 transport 图的外部独立验证。
print('=== PF smoothness on joint atlas（同源一致性检查） ===', flush=True)
pf = pd.read_parquet(root/'features/transport/n_transport_features_v1.parquet')[['jid','PF_mean']]
eps = 1e-3
for carrier, cfg, pf_file in [('n',['Eg','m_electron','n_transport'],'n_transport_features_v1.parquet'),
                              ('p',['Eg','m_hole','p_transport'],'p_transport_features_v1.parquet')]:
    layers = layers_for(cfg)
    r = embed_and_eval(layers, LAM, False)
    cons = r['coords'][:len(all_jids)]
    pff = pd.read_parquet(root/'features/transport'/pf_file)[['jid','PF_mean']].set_index('jid')
    # 只取有 transport 的 consensus
    tj = jids_map['n_transport_v1'] if carrier=='n' else jids_map['p_transport_v1']
    idx = [all_jids.index(j) for j in tj if j in all_jids]
    C = cons[idx]
    y = np.log10(np.maximum(pff.reindex(tj)['PF_mean'].values, eps))
    # 用 consensus kNN 图
    Dc = np.linalg.norm(C[:,None,:]-C[None,:,:], axis=-1)
    order = np.argsort(Dc, axis=1)[:,1:16]
    sm = np.mean([ (y[i]-y[j])**2 for i in range(len(C)) for j in order[i] ])
    # null
    null = np.array([np.mean([(y[i]-y[rng.permutation(len(C))[0]])**2 for i in range(len(C))]) for _ in range(200)])
    z = (sm-null.mean())/(null.std()+1e-12)
    p = (null<=sm).mean()
    print(f"  {carrier} joint PF smoothness={sm:.4f} z={z:.1f} p={p:.4f}（同源一致性，非独立验证）", flush=True)
    # 保存 consensus coords
    np.save(root/'data/processed'/(f'consensus_{carrier}_property.npy'), cons)

print('done', flush=True)
