"""高 ZT_e 流形 + 参数关联验证（在 Phase P-S 的 consensus 空间上）。"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))

root = Path(__file__).resolve().parents[1]
L = 2.44e-8
formula = pd.read_parquet(root/'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
edf = pd.read_parquet(root/'features/electronic/electronic_features_v1.parquet').set_index('jid')

for carrier in ['n','p']:
    cons = pd.read_parquet(root/'manifolds'/f'{carrier}_atlas_consensus.parquet')
    tdf = pd.read_parquet(root/f'features/transport/{carrier}_transport_features_v1.parquet').set_index('jid')
    df = cons.join(tdf[['S_median','log_sigma_dom_geo','PF_mean']], on='jid')
    df = df.join(edf[['Eg_optb88vdw','m_elec_median','m_hole_median']], on='jid')
    S_V = df['S_median'].values * 1e-6
    df['ZT_e'] = S_V**2 / L
    df['logPF'] = np.log10(np.maximum(df['PF_mean'], 1e-3))
    df['absS'] = df['S_median'].abs()
    df['metal'] = (df['Eg_optb88vdw'] == 0).astype(int)
    print(f'\n=== {carrier}-type consensus manifold 上的参数关联 ===')
    # Φ1, Φ2 与参数的相关
    for ph in ['Phi_1','Phi_2','Phi_3']:
        line = f'  {ph}:'
        for x in ['ZT_e','logPF','absS','log_sigma_dom_geo','Eg_optb88vdw','m_elec_median','m_hole_median']:
            sub = df[[ph,x]].dropna()
            if len(sub) > 30:
                r = stats.spearmanr(sub[ph], sub[x])[0]
                line += f' {x}={r:+.2f}'
        print(line)
    # ZT_e 高/低 在 Φ1-Φ2 空间的质心
    hi = df[df['ZT_e'] >= df['ZT_e'].quantile(0.8)]
    lo = df[df['ZT_e'] <= df['ZT_e'].quantile(0.2)]
    print(f'  high-ZT_e region centroid: Phi1={hi.Phi_1.mean():.3f} Phi2={hi.Phi_2.mean():.3f} Phi3={hi.Phi_3.mean():.3f}')
    print(f'  low-ZT_e  region centroid: Phi1={lo.Phi_1.mean():.3f} Phi2={lo.Phi_2.mean():.3f} Phi3={lo.Phi_3.mean():.3f}')
    print(f'  high vs low Eg: {hi.Eg_optb88vdw.median():.2f} vs {lo.Eg_optb88vdw.median():.2f} | m_elec: {hi.m_elec_median.median():.2f} vs {lo.m_elec_median.median():.2f}')
    df.to_parquet(root/f'manifolds/{carrier}_atlas_consensus_zte.parquet', index=False)

# 保存合并 n+p ZT_e 表
print('\n=== 保存 ZT_e 总表 ===')
allrows = []
for carrier in ['n','p']:
    tdf = pd.read_parquet(root/f'features/transport/{carrier}_transport_features_v1.parquet').set_index('jid')
    df = tdf.join(edf[['Eg_optb88vdw','m_elec_median','m_hole_median']], how='inner')
    S_V = df['S_median'].values*1e-6
    df['ZT_e'] = S_V**2/L
    df = df.reset_index()
    df['carrier'] = carrier
    df['formula'] = df['jid'].map(formula)
    allrows.append(df[['jid','carrier','formula','ZT_e','S_median','log_sigma_dom_geo','PF_mean','Eg_optb88vdw','m_elec_median','m_hole_median']])
zt = pd.concat(allrows)
zt.to_csv(root/'data/processed/ZT_e_all.csv', index=False)
print(f'ZT_e_all.csv: {zt.shape}, ZT_e>1: {(zt.ZT_e>1).sum()}, ZT_e>2: {(zt.ZT_e>2).sum()}')
