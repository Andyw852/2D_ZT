"""设计规律: 现有 JARVIS dft_2d 的 PF 决定因素 + high-PF 区域定位。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))

root = Path(__file__).resolve().parents[1]
formula = pd.read_parquet(root/'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
edf = pd.read_parquet(root/'features/electronic/electronic_features_v1.parquet').set_index('jid')

for carrier in ['n','p']:
    tdf = pd.read_parquet(root/f'features/transport/{carrier}_transport_features_v1.parquet').set_index('jid')
    df = tdf.join(edf[['Eg_optb88vdw','m_elec_median','m_hole_median']], how='inner')
    df['logPF'] = np.log10(np.maximum(df['PF_mean'], 1e-3))
    df['logsigma'] = df['log_sigma_dom_geo']
    df['absS'] = df['S_median'].abs()
    print(f'\n=== {carrier}-type (N={len(df)}) PF 决定因素 ===')
    for x in ['absS','logsigma','Eg_optb88vdw','m_elec_median','m_hole_median','S_median','S_MAD','D_sigma','A_sigma_dom']:
        sub = df[[x,'logPF']].dropna()
        rp = stats.pearsonr(sub[x], sub['logPF'])[0]
        rs = stats.spearmanr(sub[x], sub['logPF'])[0]
        print(f'  logPF vs {x:<16}: pearson={rp:+.3f}  spearman={rs:+.3f}')

    # 高 PF (top 20%) vs 低 PF (bottom 20%) 参数对比
    q80 = df['logPF'].quantile(0.8); q20 = df['logPF'].quantile(0.2)
    hi = df[df['logPF'] >= q80]; lo = df[df['logPF'] <= q20]
    print(f'\n  --- high-PF (top20%) vs low-PF (bottom20%) ---')
    for x in ['absS','logsigma','Eg_optb88vdw','m_elec_median','D_sigma']:
        a = hi[x].median(); b = lo[x].median()
        print(f'    {x:<16}: high={a:.3f}  low={b:.3f}  ratio_high/low={a/(b+1e-9):.2f}')

    # PF 峰值区域: 按 absS 与 logsigma 二维分桶
    print(f'\n  --- PF 在 (absS, logsigma) 空间的峰值 (中位 logPF) ---')
    df['S_bin'] = pd.qcut(df['absS'], 5, labels=['S1','S2','S3','S4','S5'])
    df['g_bin'] = pd.qcut(df['logsigma'], 5, labels=['g1','g2','g3','g4','g5'])
    piv = df.pivot_table(index='S_bin', columns='g_bin', values='logPF', aggfunc='median')
    print(piv.round(2).to_string())
    # 最高 PF 材料 top 10
    top = df.sort_values('PF_mean', ascending=False).head(10)
    print(f'\n  --- top-10 PF 材料 ---')
    for j, r in top.iterrows():
        print(f'    {j} ({formula.get(j,"")}): PF={r["PF_mean"]:.1f} S={r["S_median"]:.1f} logsigma={r["logsigma"]:.2f} Eg={r["Eg_optb88vdw"]:.2f}')

# 保存 n/p 合并表（供后续）
for carrier in ['n','p']:
    tdf = pd.read_parquet(root/f'features/transport/{carrier}_transport_features_v1.parquet').set_index('jid')
    df = tdf.join(edf, how='inner')
    df.to_csv(root/f'data/processed/PF_analysis_{carrier}.csv', index=True)
print('\nsaved PF_analysis_{n,p}.csv')
