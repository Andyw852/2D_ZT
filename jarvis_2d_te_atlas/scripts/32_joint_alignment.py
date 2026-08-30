"""Phase P-Q: joint alignment —— λ 扫描 + structure inclusion + random anchor + layer ablation。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from multiview_utils import build_supra, joint_embedding, load_W, jids_of, scale_layer
from scipy import sparse

root = Path(__file__).resolve().parents[1]
gdir = root / 'graphs'
K = 15
N_EIGS = 20
rng = np.random.RandomState(42)

# 载入 scaled layers
d = np.load(gdir / 'multiview_scaled.npz', allow_pickle=True)
scaled = {k: d[k].item() if hasattr(d[k],'item') else d[k] for k in d.files}
jids_map = {name: jids_of(name, gdir) for name in ['structure_v1','Eg_v1','m_electron_v1','m_hole_v1','n_transport_v1','p_transport_v1']}
# 修正 mass jids 从 mass graph nodes csv
jids_map['m_electron_v1'] = jids_of('m_electron_v1', gdir)
jids_map['m_hole_v1'] = jids_of('m_hole_v1', gdir)

all_jids = sorted(jids_map['structure_v1'])

def make_layers(config):
    mapping = {'structure':'structure_v1','Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1',
               'n_transport':'n_transport_v1','p_transport':'p_transport_v1'}
    return [(mapping[v], scaled[mapping[v]], jids_map[mapping[v]]) for v in config]

def eval_config(layers, lam, consensus_jids, n_eigs=N_EIGS):
    A, offsets, cid = build_supra(layers, consensus_jids, lam)
    ev, V = joint_embedding(A, n_eigs)
    if V is None:
        return None
    coords = V[:, 1:n_eigs+1]  # 去掉 trivial
    # consensus coords
    Nc = len(consensus_jids)
    cons = coords[:Nc]
    # tension（同 JID 的 consensus 与 view copy 距离，归一化）
    all_t = []
    for name, W, jids in layers:
        jmap = {j:i for i,j in enumerate(jids)}
        off = offsets[name]
        for jid in jids:
            if jid in cid:
                ci = cid[jid]; vi = off + jmap[jid]
                all_t.append(np.linalg.norm(cons[ci] - coords[vi]))
    all_t = np.array(all_t)
    # 归一化：除以所有 joint pair 距离的中位数（采样）
    n_tot = coords.shape[0]
    sample = min(n_tot, 3000)
    idx = rng.choice(n_tot, sample, replace=False)
    cc = coords[idx]
    pair_d = np.linalg.norm(cc[:,None,:] - cc[None,:,:], axis=-1)
    median_pair = np.median(pair_d[np.triu_indices(sample, k=1)])
    tension_norm = all_t / (median_pair + 1e-12)
    # preservation
    pres = {}
    for name, W, jids in layers:
        # original kNN from W
        jmap = {j:i for i,j in enumerate(jids)}
        Wc = W.tocsr()
        orig = {}
        for i in range(W.shape[0]):
            row = Wc[i].toarray().ravel()
            orig[jids[i]] = set(np.argsort(-row)[:K+1].tolist()) - {i}
        # joint kNN from view coords
        off = offsets[name]
        vc = coords[off:off+len(jids)]
        Dv = np.linalg.norm(vc[:,None,:] - vc[None,:,:], axis=-1)
        order = np.argsort(Dv, axis=1)[:, 1:K+1]
        joint = {jids[i]: set(order[i].tolist()) for i in range(len(jids))}
        ov = np.mean([len(orig[j] & joint[j]) / K for j in jids])
        pres[name] = ov
    # giant component / spectral gap
    ncomp, lab = sparse.csgraph.connected_components(A, directed=False)
    sz = np.bincount(lab)
    giant = sz.max() / A.shape[0]
    spectral_gap = float(ev[1]) if len(ev) > 1 else np.nan  # 第一非平凡本征值
    return {'median_tension': float(np.median(tension_norm)), 'mean_tension': float(tension_norm.mean()),
            'P': pres, 'giant': giant, 'spectral_gap': spectral_gap, 'n_nodes': A.shape[0],
            'coords': coords, 'offsets': offsets, 'cid': cid}

configs = {
    'n_Full': ['structure','Eg','m_electron','n_transport'],
    'n_Property': ['Eg','m_electron','n_transport'],
    'p_Full': ['structure','Eg','m_hole','p_transport'],
    'p_Property': ['Eg','m_hole','p_transport'],
}
LAMBDAS = [0.01, 0.03, 0.10, 0.30, 1.00, 3.00, 10.00]

scan_rows = []
coords_cache = {}
for cname, cfg in configs.items():
    layers = make_layers(cfg)
    carrier = 'n' if 'n_transport' in cfg else 'p'
    prev_cons = None
    for lam in LAMBDAS:
        r = eval_config(layers, lam, all_jids)
        if r is None:
            continue
        row = {'config': cname, 'carrier': carrier, 'lambda': lam, 'median_anchor_tension': round(r['median_tension'],4),
               'giant_component_fraction': round(r['giant'],4), 'spectral_gap': round(r['spectral_gap'],4)}
        for view in ['structure','Eg','m_electron','m_hole','n_transport','p_transport']:
            key = {'structure':'structure_v1','Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1','n_transport':'n_transport_v1','p_transport':'p_transport_v1'}.get(view)
            row['P_' + view] = round(r['P'].get(key, np.nan), 4)
        # consensus stability vs previous lambda
        cons = r['coords'][:len(all_jids)]
        if prev_cons is not None:
            D1 = np.linalg.norm(cons[:,None,:]-cons[None,:,:], axis=-1); D2 = np.linalg.norm(prev_cons[:,None,:]-prev_cons[None,:,:], axis=-1)
            o1 = np.argsort(D1, axis=1)[:,1:16]; o2 = np.argsort(D2, axis=1)[:,1:16]
            stab = np.mean([len(set(o1[i])&set(o2[i]))/15 for i in range(len(cons))])
            row['consensus_stability'] = round(stab, 4)
        else:
            row['consensus_stability'] = np.nan
        scan_rows.append(row)
        prev_cons = cons
        coords_cache[(cname, lam)] = (r['coords'], r['offsets'], r['cid'])
        print(f"{cname} lambda={lam}: tension={row['median_anchor_tension']:.3f} P_transport={row.get('P_n_transport',row.get('P_p_transport',np.nan))} P_structure={row.get('P_structure',np.nan)} P_Eg={row.get('P_Eg',np.nan)} P_mass={row.get('P_m_electron',row.get('P_m_hole',np.nan))}")

sdf = pd.DataFrame(scan_rows)
sdf.to_csv(root / 'data/audit/lambda_scan_all.csv', index=False)
for carrier in ['n','p']:
    sdf[sdf['carrier']==carrier].to_csv(root / 'data/audit' / ('lambda_scan_' + carrier + '.csv'), index=False)
print('\nsaved lambda_scan_n.csv / lambda_scan_p.csv')

# 保存 coords cache 供后续（存为 npz）
np.savez_compressed(root / 'data/processed/joint_coords_cache.npz', 
                    **{('_'.join(map(str,k))): v[0] for k,v in coords_cache.items()}, allow_pickle=True)
print('saved joint_coords_cache.npz')
