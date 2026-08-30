import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
sys.path.insert(0, 'scripts')
from multiview_utils import build_supra, joint_embedding, load_W, jids_of

root = Path('.')
gdir = root / 'graphs'
N_EIGS = 20; LAM = 0.3
rng = np.random.RandomState(42)
d = np.load(gdir / 'multiview_scaled.npz', allow_pickle=True)
scaled = {k: d[k].item() if hasattr(d[k], 'item') else d[k] for k in d.files}
jids_map = {n: jids_of(n, gdir) for n in ['structure_v1','Eg_v1','m_electron_v1','m_hole_v1','n_transport_v1','p_transport_v1']}
all_jids = sorted(jids_map['structure_v1'])

def layers_for(cfg):
    m = {'Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1','n_transport':'n_transport_v1','p_transport':'p_transport_v1'}
    return [(m[v], scaled[m[v]], jids_map[m[v]]) for v in cfg]

def abs_tension(layers, lam, shuffle=False):
    layers2 = layers
    if shuffle:
        layers2 = []
        for name, W, jids in layers:
            perm = rng.permutation(len(jids)); jids2 = [jids[p] for p in perm]
            layers2.append((name, W, jids2))
    A, offsets, cid = build_supra(layers2, all_jids, lam)
    ev, V = joint_embedding(A, N_EIGS)
    coords = V[:, 1:N_EIGS+1]
    cons = coords[:len(all_jids)]
    ts = []
    for name, W, jids in layers2:
        jmap = {j: i for i, j in enumerate(jids)}
        off = offsets[name]
        for jid in jids:
            if jid in cid:
                ts.append(np.linalg.norm(cons[cid[jid]] - coords[off + jmap[jid]]))
    return float(np.median(np.array(ts)))

rows = []
for carrier, cfg in [('n', ['Eg','m_electron','n_transport']), ('p', ['Eg','m_hole','p_transport'])]:
    layers = layers_for(cfg)
    real = abs_tension(layers, LAM, shuffle=False)
    null = np.array([abs_tension(layers, LAM, shuffle=True) for _ in range(60)])
    z = (real - null.mean()) / (null.std() + 1e-12)
    p = (null <= real).mean()
    rows.append({'carrier': carrier, 'real_tension': round(real,5), 'null_mean': round(float(null.mean()),5),
                 'null_std': round(float(null.std()),5), 'z': round(float(z),2), 'p': round(float(p),4)})
    print(f"{carrier}: real={real:.5f} null={null.mean():.5f}+-{null.std():.5f} z={z:.1f} p={p:.4f}", flush=True)
pd.DataFrame(rows).to_csv(root / 'data/audit/random_anchor_abs.csv', index=False)
print('saved random_anchor_abs.csv', flush=True)
