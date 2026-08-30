"""晶格参数预筛配对（非超晶格构造）：晶格失配 + 输运对比。

严格说本脚本只做原胞 a、b、area、晶格夹角 ang 的预筛，不生成双层 POSCAR、
不搜索共格超胞/旋转/交换晶轴、不考虑扭角/层间距/堆垛/界面终止/应变能。
输出应理解为"晶格参数预筛配对"，而非"已构造的超晶格候选结构"。
"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

root = Path(__file__).resolve().parents[1]
sdf = pd.read_parquet(root/'data/processed/standardized_2d_structures.parquet').set_index('jid')
formula = sdf['formula'].to_dict()
zt = pd.read_csv(root/'data/processed/ZT_e_all.csv')

def lattice_params(jid):
    cell = np.array(json.loads(sdf.loc[jid,'lattice']))
    a = np.linalg.norm(cell[0]); b = np.linalg.norm(cell[1])
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    ang = np.degrees(np.arccos(np.clip(np.dot(cell[0],cell[1])/(a*b+1e-12), -1, 1)))
    return a, b, area, ang

# 全部材料的晶格参数
params = {j: lattice_params(j) for j in sdf.index}
# 晶格失配: a/b、面积、夹角
def mismatch(j1, j2):
    a1,b1,area1,ang1 = params[j1]; a2,b2,area2,ang2 = params[j2]
    ma = abs(a1-a2)/min(a1,a2); mb = abs(b1-b2)/min(b1,b2)
    marea = abs(area1-area2)/min(area1,area2)
    mang = abs(ang1-ang2)  # 晶格夹角差(度)
    return max(ma, mb), marea, ma, mb, mang

# 每 carrier 取 top-40 ZT_e 材料
THRESH_LATTICE = 0.05  # 5% 晶格失配
THRESH_ANGLE   = 5.0   # 晶格夹角差(度)
print('=== 晶格参数预筛配对（晶格失配<5% 且 夹角差<5°）===')
pairs = []
seen = set()  # 无序去重：排除 (A,B)/(B,A) 方向重复
for carrier in ['n','p']:
    top = zt[zt['carrier']==carrier].sort_values('ZT_e', ascending=False).head(40)
    other_carrier = 'p' if carrier=='n' else 'n'
    for _, row in top.iterrows():
        jA = row['jid']
        # 找相反 carrier 的高 ZT_e 材料作为异质结伙伴
        pool = zt[zt['carrier']==other_carrier].sort_values('ZT_e', ascending=False).head(60)
        for _, r2 in pool.iterrows():
            jB = r2['jid']
            if jA == jB:      # 排除自配对（同一 JID 在 n/p 两池都出现）
                continue
            key = tuple(sorted((str(jA), str(jB))))
            if key in seen:   # 去方向重复
                continue
            seen.add(key)
            m_lat, m_area, ma, mb, mang = mismatch(jA, jB)
            if m_lat < THRESH_LATTICE and mang < THRESH_ANGLE:
                pairs.append({
                    'A_jid': jA, 'A_formula': formula.get(jA,''), 'A_carrier': carrier, 'A_ZT_e': round(row['ZT_e'],2),
                    'A_S': round(row['S_median'],1), 'A_Eg': round(row['Eg_optb88vdw'],2),
                    'B_jid': jB, 'B_formula': formula.get(jB,''), 'B_carrier': other_carrier, 'B_ZT_e': round(r2['ZT_e'],2),
                    'B_S': round(r2['S_median'],1), 'B_Eg': round(r2['Eg_optb88vdw'],2),
                    'lattice_mismatch_pct': round(m_lat*100,2), 'area_mismatch_pct': round(m_area*100,2),
                    'angle_mismatch_deg': round(mang,1),
                    'S_contrast': round(abs(row['S_median']-r2['S_median']),1),
                })
pdf = pd.DataFrame(pairs)
pdf = pdf.sort_values(['lattice_mismatch_pct','A_ZT_e'], ascending=[True, False])
pdf.to_csv(root/'data/processed/superlattice_candidate_pairs.csv', index=False)
print(f'晶格参数预筛配对（无序、排除自配对）: {len(pdf)} 对')
print(pdf.head(25).to_string(index=False))
print(f'\n保存 superlattice_candidate_pairs.csv（预筛配对，{len(pdf)} 对；非已构造超晶格）')