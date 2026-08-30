"""电子 ZT 上限与高-ZT 天花板定位（含单位核对与定义性说明）。

ZT_e^WF = S^2/L（固定金属 Lorenz 常数 L=2.44e-8 W·Ω/K²）是「上限假设」：
它只含 Seebeck S，因此「高 ZT_e 由高 S 主导」是定义使然（循环），不是独立发现。
真实电子 ZT = S^2·σ·T/κ_e 需要 σ 与 κ_e 数据；但 JARVIS 的 nkappa 与 ncond 缩放
不一致（见下方单位核对），不能直接相除。最终真实 ZT 还需声子 κ_L（MACE 计算）。"""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

root = Path(__file__).resolve().parents[1]
L = 2.44e-8  # 金属极限 Lorenz 常数 W·Ω/K^2（高 Seebeck/半导体下实际 L 会偏离）
formula = pd.read_parquet(root/'data/processed/standardized_2d_structures.parquet').set_index('jid')['formula'].to_dict()
edf = pd.read_parquet(root/'features/electronic/electronic_features_v1.parquet').set_index('jid')

for carrier in ['n','p']:
    tdf = pd.read_parquet(root/f'features/transport/{carrier}_transport_features_v1.parquet').set_index('jid')
    df = tdf.join(edf[['Eg_optb88vdw','m_elec_median','m_hole_median']], how='inner')
    # S in V/K (S_median 是 μV/K)
    S_V = df['S_median'].values * 1e-6
    # 电子 ZT 上限（固定金属 L 的 WF 上限）：ZT_e^WF = S^2/L。
    # 注意：此定义下 ZT_e 只含 S，故「ZT_e 与 S 相关」是定义使然，非独立验证。
    ZTe = S_V**2 / L
    df['ZT_e_ceiling'] = ZTe
    # 【单位核对】ZT_e = S^2·σ·T/κ_e 需要 σ、κ_e 原始值。经 PF 自洽验证：
    #   ncond 单位 S/m（PF_raw = S_μV²·σ_raw×1e-6 自洽成立）；
    #   但 nkappa 与 ncond 缩放不一致：κ_raw/(L·σ_raw·T) ≈ 1e14（非恒定，材料间 9e13~5e14），
    #   即 JARVIS 的 nkappa 不能直接与 ncond 相除算 ZT_e，需先校准缩放。
    # 因此在未校准 nkappa 前，只用 ZT_e^WF = S^2/L 作「上限」参考，并明确其定义性。
    T = 600.0
    print(f'\n=== {carrier}-type 电子 ZT 上限 ZT_e = S^2/L ===')
    print(f'  ZT_e 分布: median={df["ZT_e_ceiling"].median():.2f} p90={df["ZT_e_ceiling"].quantile(0.9):.2f} max={df["ZT_e_ceiling"].max():.2f}')
    print(f'  ZT_e > 1 的材料数: {(df["ZT_e_ceiling"]>1).sum()} / {len(df)}')
    print(f'  ZT_e > 2 的材料数: {(df["ZT_e_ceiling"]>2).sum()} / {len(df)}')
    # 高 ZT_e 材料的 Eg / S / m* 特征
    top = df.sort_values('ZT_e_ceiling', ascending=False).head(12)
    print(f'\n  --- 高 ZT_e 天花板 top-12（注：真实 ZT = ZT_e * kappa_e/(kappa_e+kappa_L) < ZT_e）---')
    print(f'  {"jid":<14}{"formula":<10}{"ZT_e":<7}{"S_median":<10}{"Eg":<7}{"m_median"}')
    for j, r in top.iterrows():
        m = r['m_elec_median'] if carrier=='n' else r['m_hole_median']
        print(f'  {j:<14}{formula.get(j,"")[:9]:<10}{r["ZT_e_ceiling"]:<7.2f}{r["S_median"]:<10.1f}{r["Eg_optb88vdw"]:<7.2f}{m:.3f}')
    # ZT_e 与各参数的相关
    from scipy import stats
    for x in ['Eg_optb88vdw','m_elec_median','m_hole_median','S_MAD','D_sigma']:
        sub = df[[x,'ZT_e_ceiling']].dropna()
        print(f'  ZT_e vs {x:<15}: spearman={stats.spearmanr(sub[x], sub["ZT_e_ceiling"])[0]:+.3f}')
    print('  （注：ZT_e^WF ∝ S²，与 |S|/S_MAD 的相关是定义使然；与 Eg/m* 的相关才含物理信息）')

# 结论性统计
n = pd.read_csv(root/'data/processed/PF_analysis_n.csv').set_index('jid')
p = pd.read_csv(root/'data/processed/PF_analysis_p.csv').set_index('jid')
print('\n=== ZT_e 上限 vs PF 的关系（定义性说明） ===')
print('（ZT_e^WF ∝ S²，PF ∝ S²σ：高 ZT_e^WF 仅由高 |S| 决定，是定义使然；')
print('  高 PF 还需高 σ，才是数据自洽的量。真实 ZT = S²σT/(κ_e+κ_L)，')
print('  其中 κ_L 由 MACE 声子计算给出，是独立于 JARVIS 输运标签的验证。）')
