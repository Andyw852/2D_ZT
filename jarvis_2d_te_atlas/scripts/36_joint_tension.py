"""Phase R-S: 最终 joint tension + consensus parquet + metal/semiconductor。"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
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
edf = pd.read_parquet(root / 'features/electronic/electronic_features_v1.parquet').set_index('jid')

def layers_for(cfg):
    m = {'Eg':'Eg_v1','m_electron':'m_electron_v1','m_hole':'m_hole_v1','n_transport':'n_transport_v1','p_transport':'p_transport_v1'}
    return [(m[v], scaled[m[v]], jids_map[m[v]]) for v in cfg]

def full_joint(cfg, carrier):
    layers = layers_for(cfg)
    A, offsets, cid = build_supra(layers, all_jids, LAM)
    ev, V = joint_embedding(A, N_EIGS)
    coords = V[:, 1:N_EIGS+1]
    cons = coords[:len(all_jids)]
    # tension per-JID per-view
    rows = []
    for jid in all_jids:
        row = {'jid': jid, 'formula': edf.reindex([jid]).index[0] if jid in edf.index else jid}
        # formula from standardized
        fs = None
        for name, W, jids in layers:
            key = name
            jmap = {j:i for i,j in enumerate(jids)}
            off = offsets[name]
            if jid in jmap:
                T = float(np.linalg.norm(cons[cid[jid]] - coords[off + jmap[jid]]))
                row['T_' + {'Eg_v1':'Eg','m_electron_v1':'mass','m_hole_v1':'mass','n_transport_v1':'transport','p_transport_v1':'transport'}[name]] = T
        tvals = [v for k,v in row.items() if k.startswith('T_')]
        if tvals:
            row['T_mean'] = float(np.mean(tvals)); row['T_max'] = float(np.max(tvals))
        row['view_count'] = len(tvals)
        rows.append(row)
    tdf = pd.DataFrame(rows)
    # metal/semiconductor
    tdf['Eg'] = tdf['jid'].map(lambda j: edf.reindex([j])['Eg_optb88vdw'].values[0] if j in edf.index else np.nan)
    tdf['metal'] = tdf['Eg'] == 0
    return tdf, cons, offsets, cid, layers

formula = pd.read_parquet(root/'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
for carrier, cfg, pf_file in [('n',['Eg','m_electron','n_transport'],'n_transport_features_v1.parquet'),
                              ('p',['Eg','m_hole','p_transport'],'p_transport_features_v1.parquet')]:
    tdf, cons, offsets, cid, layers = full_joint(cfg, carrier)
    tdf['formula'] = tdf['jid'].map(formula)
    tdf.to_csv(root / 'data/processed' / ('joint_tension_' + carrier + '.csv'), index=False)
    # consensus parquet
    cons_df = pd.DataFrame(cons, columns=[f'Phi_{i}' for i in range(1, N_EIGS+1)])
    cons_df.insert(0, 'jid', all_jids)
    # add meta
    transport_jids = set(jids_map['n_transport_v1'] if carrier=='n' else jids_map['p_transport_v1'])
    mass_jids = set(jids_map['m_electron_v1'] if carrier=='n' else jids_map['m_hole_v1'])
    cons_df['view_count'] = 1 + cons_df['jid'].isin(set(jids_map['Eg_v1'])).astype(int) + cons_df['jid'].isin(mass_jids).astype(int) + cons_df['jid'].isin(transport_jids).astype(int)
    cons_df['transport_informed'] = cons_df['jid'].isin(transport_jids)
    cons_df['mass_available'] = cons_df['jid'].isin(mass_jids)
    cons_df.to_parquet(root / 'manifolds' / (carrier + '_atlas_consensus.parquet'), index=False)
    # metal/semiconductor 统计
    metals = tdf[tdf['metal']==True]; semis = tdf[tdf['metal']==False]
    print(f'{carrier}: tension rows={len(tdf)}, metals={len(metals)}, semis={len(semis)}')
    print(f'  T_transport: median={tdf["T_transport"].median():.4f} (n={tdf["T_transport"].notna().sum()})')
    print(f'  T_Eg: median={tdf["T_Eg"].median():.4f}')
    print(f'  top-5 T_transport JIDs:', tdf.sort_values('T_transport', ascending=False).head(5)['jid'].tolist())
    print(f'  top-5 T_Eg JIDs:', tdf.sort_values('T_Eg', ascending=False).head(5)['jid'].tolist())
print('saved joint_tension + consensus parquet')
