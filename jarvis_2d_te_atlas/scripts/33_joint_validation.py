"""Phase Q-R-S: random anchor + layer ablation + duplicate sensitivity + tension + PF + n/p 比较。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
sys.path.insert(0, str(Path(__file__).resolve().parent))
from multiview_utils import build_supra, joint_embedding, load_W, jids_of, scale_layer

root = Path(__file__).resolve().parents[1]
gdir = root / 'graphs'
K = 15
N_EIGS = 20
LAM = 0.3
rng = np.random.RandomState(42)

d = np.load(gdir / 'multiview_scaled.npz', allow_pickle=True)
scaled = {k: d[k].item() if hasattr(d[k],'item') else d[k] for k in d.files}
jids_map = {n: jids_of(n, gdir) for n in ['structure_v1','Eg_v1','m_electron_v1','m_hole_v1','n_transport_v1','p_transport_v1']}
all_jids = sorted(jids_map['structure_v1'])

def layers_for(config):
    m = {'structure':'structure_v1','Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1','n_transport':'n_transport_v1','p_transport':'p_transport_v1'}
    return [(m[v], scaled[m[v]], jids_map[m[v]]) for v in config]

def eval_joint(layers, lam, consensus_jids, n_eigs=N_EIGS, shuffle=False):
    A, offsets, cid = build_supra(layers, consensus_jids, lam)
    if shuffle:
        # 打乱 identity edges: 重排每个 view 的 jid 顺序
        layers2 = []
        for name, W, jids in layers:
            perm = rng.permutation(len(jids))
            jids2 = [jids[p] for p in perm]
            # W 的行对应原 jids 顺序，identity 用打乱后的 jids2 但 W 不变
            layers2.append((name, W, jids2))
        A, offsets, cid = build_supra(layers2, consensus_jids, lam)
    ev, V = joint_embedding(A, n_eigs)
    if V is None:
        return None
    coords = V[:, 1:n_eigs+1]
    Nc = len(consensus_jids)
    cons = coords[:Nc]
    # tension
    ts = []
    for name, W, jids in (layers2 if shuffle else layers):
        jmap = {j:i for i,j in enumerate(jids)}
        off = offsets[name]
        for jid in jids:
            if jid in cid:
                ts.append(np.linalg.norm(cons[cid[jid]] - coords[off+jmap[jid]]))
    ts = np.array(ts)
    # normalized
    n_tot = coords.shape[0]; sample = min(n_tot, 3000)
    idx = rng.choice(n_tot, sample, replace=False)
    cc = coords[idx]; pd_ = np.linalg.norm(cc[:,None,:]-cc[None,:,:], axis=-1)
    mp = np.median(pd_[np.triu_indices(sample, k=1)])
    return {'median_tension': float(np.median(ts/(mp+1e-12))), 'coords': coords, 'offsets': offsets, 'cid': cid,
            'layers': (layers2 if shuffle else layers)}

# 1) random anchor control (n/p property, lam=0.3)
print('=== random anchor control (lam=0.3, 100 perms) ===')
ra_rows = []
for carrier, cfg in [('n', ['Eg','m_electron','n_transport']), ('p', ['Eg','m_hole','p_transport'])]:
    layers = layers_for(cfg)
    real = eval_joint(layers, LAM, all_jids, shuffle=False)
    null = [eval_joint(layers, LAM, all_jids, shuffle=True)['median_tension'] for _ in range(100)]
    null = np.array(null)
    z = (real['median_tension'] - null.mean()) / (null.std() + 1e-12)
    p = (null <= real['median_tension']).mean()
    ra_rows.append({'carrier': carrier, 'real_tension': round(real['median_tension'],4),
                    'null_mean': round(null.mean(),4), 'null_std': round(null.std(),4),
                    'z_score': round(z,2), 'empirical_p': round(p,4)})
    print(f'  {carrier}: real={real["median_tension"]:.4f} null={null.mean():.4f}+-{null.std():.4f} z={z:.1f} p={p:.4f}')
pd.DataFrame(ra_rows).to_csv(root / 'data/audit/random_anchor.csv', index=False)

# 2) layer ablation (n/p property)
print('=== layer ablation (lam=0.3) ===')
ab_rows = []
for carrier, cfg in [('n', ['Eg','m_electron','n_transport']), ('p', ['Eg','m_hole','p_transport'])]:
    full = eval_joint(layers_for(cfg), LAM, all_jids)
    full_cons = full['coords'][:len(all_jids)]
    for drop in cfg:
        cfg2 = [v for v in cfg if v != drop]
        r = eval_joint(layers_for(cfg2), LAM, all_jids)
        cons = r['coords'][:len(all_jids)]
        # consensus kNN overlap between full and ablated
        D1 = np.linalg.norm(full_cons[:,None,:]-full_cons[None,:,:], axis=-1)
        D2 = np.linalg.norm(cons[:,None,:]-cons[None,:,:], axis=-1)
        o1 = np.argsort(D1,axis=1)[:,1:16]; o2 = np.argsort(D2,axis=1)[:,1:16]
        ov = np.mean([len(set(o1[i])&set(o2[i]))/15 for i in range(len(cons))])
        ab_rows.append({'carrier': carrier, 'dropped': drop, 'consensus_knn_overlap': round(ov,4), 'median_tension': round(r['median_tension'],4)})
        print(f'  {carrier} drop {drop}: consensus overlap={ov:.4f} tension={r["median_tension"]:.4f}')
pd.DataFrame(ab_rows).to_csv(root / 'data/audit/layer_ablation.csv', index=False)

# 3) duplicate sensitivity: 用 collapsed 代表 JID 真正重跑联合流形，比较 transport 邻域保持
print('=== duplicate sensitivity (collapsed re-run) ===')
dstruct = np.load(root / 'data/processed/d_struct_baseline.npy')
sjids = jids_map['structure_v1']
dup_groups = {}
iu = np.triu_indices(len(sjids), k=1)
for t in range(len(iu[0])):
    if dstruct[iu[0][t], iu[1][t]] == 0:
        a, b = sjids[iu[0][t]], sjids[iu[1][t]]
        ga, gb = dup_groups.get(a, {a}), dup_groups.get(b, {b})
        g = ga | gb
        for x in g: dup_groups[x] = g
groups = set(map(tuple, dup_groups.values()))
n_dup = sum(1 for g in groups if len(g) > 1)
rep = {}
for g in groups:
    for x in g: rep[x] = sorted(g)[0]
collapsed_jids = sorted({rep.get(j, j) for j in all_jids})  # 重复组内只留代表，非重复材料保留自身
print(f'  exact duplicate groups: {n_dup}, collapsed consensus JIDs: {len(collapsed_jids)} (vs full {len(all_jids)})')

# 真正重跑：full(全 1103) vs collapsed(去重代表) 联合流形，比较代表 JID 的 transport kNN 邻域 overlap
dup_rows = []
for carrier, cfg, tname in [('n', ['Eg','m_electron','n_transport'], 'n_transport_v1'),
                            ('p', ['Eg','m_hole','p_transport'], 'p_transport_v1')]:
    full = eval_joint(layers_for(cfg), LAM, all_jids)
    coll = eval_joint(layers_for(cfg), LAM, collapsed_jids)
    tj = jids_map[tname]
    common = [j for j in collapsed_jids if j in all_jids and j in tj]
    fi = [all_jids.index(j) for j in common]
    ci = [collapsed_jids.index(j) for j in common]
    Cfull = full['coords'][fi]
    Ccoll = coll['coords'][ci]
    Df = np.linalg.norm(Cfull[:,None,:]-Cfull[None,:,:], axis=-1)
    Dc = np.linalg.norm(Ccoll[:,None,:]-Ccoll[None,:,:], axis=-1)
    of = np.argsort(Df,axis=1)[:,1:16]; oc = np.argsort(Dc,axis=1)[:,1:16]
    ov = np.mean([len(set(of[i])&set(oc[i]))/15 for i in range(len(Cfull))])
    passed = bool(ov >= 0.9)
    dup_rows.append({'carrier': carrier, 'n_common': len(common),
                     'consensus_knn_overlap': round(ov,4), 'passed': passed})
    print(f'  {carrier}: collapsed kNN(15) overlap={ov:.4f} (N={len(common)}) passed={passed}')
pd.DataFrame(dup_rows).to_csv(root/'data/audit/duplicate_sensitivity.csv', index=False)
dup_passed = bool(len(dup_rows) and all(r['passed'] for r in dup_rows))
print(f'  DUPLICATE_SENSITIVITY_PASSED = {dup_passed}')

# 4) n/p atlas comparison（property, lam=0.3）
print('=== n/p atlas comparison ===')
rn = eval_joint(layers_for(['Eg','m_electron','n_transport']), LAM, all_jids)
rp = eval_joint(layers_for(['Eg','m_hole','p_transport']), LAM, all_jids)
common = sorted(set(jids_map['n_transport_v1']) & set(jids_map['p_transport_v1']))
cni = [all_jids.index(j) for j in common if j in all_jids]
if cni:
    Cn = rn['coords'][cni]; Cp = rp['coords'][cni]
    Dn = np.linalg.norm(Cn[:,None,:]-Cn[None,:,:], axis=-1); Dp = np.linalg.norm(Cp[:,None,:]-Cp[None,:,:], axis=-1)
    from scipy import stats
    rho = stats.spearmanr(Dn[np.triu_indices(len(cni),k=1)], Dp[np.triu_indices(len(cni),k=1)])[0]
    on = np.argsort(Dn,axis=1)[:,1:16]; op = np.argsort(Dp,axis=1)[:,1:16]
    ov = np.mean([len(set(on[i])&set(op[i]))/15 for i in range(len(cni))])
    print(f'  n vs p consensus (N={len(cni)}): distance Spearman={rho:.4f}, kNN(15) overlap={ov:.4f}')

# 保存 n/p 共同 JID 的 consensus 坐标（供后续 PF/tension）
np.save(root / 'data/processed/consensus_n_property_lam0.3.npy', rn['coords'])
np.save(root / 'data/processed/consensus_p_property_lam0.3.npy', rp['coords'])
print('saved consensus coordinates')
