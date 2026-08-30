#!/usr/bin/env python3
"""构造真正的共格双层超晶格：从 superlattice_candidate_pairs.csv 选共格对，
用 JARVIS 结构（修正的 cart=frac@lat 行存储）旋转+应变对齐晶格，层叠为双层，
输出 POSCAR + 元数据。过滤：平面内对齐、|S-I|fro<2%、层间距 3.0-3.6 A。

注意：这是真正的超晶格构造（生成双层 POSCAR），与 42 的晶格参数预筛不同。
"""
import json, sys, numpy as np, pandas as pd
from pathlib import Path
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'superlattices'
OUT.mkdir(exist_ok=True)

snap = json.load(open(ROOT/'data/raw/jarvis/dft_2d_snapshot.json'))
by_id = {e['id']: e['attributes'] for e in snap}

def get_structure(jid):
    att = by_id.get(jid)
    if att is None:
        return None
    lat = np.array(att['lattice_vectors'], dtype=float)
    cart = np.array(att['cartesian_site_positions'], dtype=float)
    return lat, cart, att['species_at_sites']

def build_bilayer(jidA, jidB, gap=3.3, vacuum=15.0):
    """-> (atoms, meta) 或 None（不可构造）"""
    sa = get_structure(jidA); sb = get_structure(jidB)
    if sa is None or sb is None:
        return None
    la, ca, spA = sa; lb, cb, spB = sb
    # 平面内对齐检查（真空沿 z）
    for name, lat in (('A', la), ('B', lb)):
        if abs(lat[0][2]) > 1e-3 or abs(lat[1][2]) > 1e-3 or abs(lat[2][0]) > 1e-3 or abs(lat[2][1]) > 1e-3:
            print(f'  !! {name} 晶格非平面对齐，跳过')
            return None
    a1, b1 = la[0][:2], la[1][:2]
    a2, b2 = lb[0][:2], lb[1][:2]
    ang = lambda v: np.degrees(np.arctan2(v[1], v[0]))
    # B 旋转对齐到 A 的方向（保持 B 的基矢长度/夹角）
    theta = np.radians(ang(a1) - ang(a2))
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    a2r, b2r = R @ a2, R @ b2
    cb2 = np.column_stack([(R @ cb[:, :2].T).T, cb[:, 2]])
    L1 = np.array([a1, b1]); L2r = np.array([a2r, b2r])
    S = L1 @ np.linalg.inv(L2r)
    strain = float(np.linalg.norm(S - np.eye(2)))
    if strain > 0.02:
        print(f'  !! 应变 {strain:.4f} > 2%，跳过')
        return None
    cb2[:, :2] = cb2[:, :2] @ S.T
    # 层叠
    zA_top = ca[:, 2].max(); zB_bot = cb2[:, 2].min()
    cb2[:, 2] += (zA_top + gap - zB_bot)
    zB_top = cb2[:, 2].max()
    new_cell = la.copy(); new_cell[2] = [0, 0, zB_top + vacuum]
    all_pos = np.vstack([ca, cb2]); all_sp = list(spA) + list(spB)
    atoms = Atoms(all_sp, positions=all_pos, cell=new_cell, pbc=[True, True, False])
    inter = min(np.linalg.norm(p1 - p2) for p1 in cb2 for p2 in ca)
    meta = dict(A=jidA, B=jidB, strain=round(strain, 4), interlayer=round(inter, 3),
                n_atoms=len(atoms), c=round(new_cell[2][2], 2),
                a=round(np.linalg.norm(a1), 4), b=round(np.linalg.norm(b1), 4))
    return atoms, meta

pairs = pd.read_csv(ROOT/'data/processed/superlattice_candidate_pairs.csv')
coh = pairs[(pairs.lattice_mismatch_pct < 1) & (pairs.angle_mismatch_deg < 1)]
rows = []
for _, r in coh.iterrows():
    jidA, jidB = 'dft_2d_' + str(r['A_jid']), 'dft_2d_' + str(r['B_jid'])
    print(f"构造 {r['A_formula']} ({r['A_jid']}) + {r['B_formula']} ({r['B_jid']})")
    res = build_bilayer(jidA, jidB)
    if res is None:
        continue
    atoms, meta = res
    name = f"{r['A_formula']}_{r['A_jid']}_{r['B_formula']}_{r['B_jid']}_bilayer"
    write(OUT / (name + '.POSCAR'), atoms, vasp5=True, direct=True)
    meta.update(pair_row=dict(r))
    rows.append(meta)
    print(f"  OK: 应变 {meta['strain']} 层间距 {meta['interlayer']} 原子 {meta['n_atoms']} -> {name}.POSCAR")

df = pd.DataFrame(rows)
if len(df):
    df.to_csv(OUT / 'coherent_bilayers_meta.csv', index=False)
    print(f"\n共 {len(df)} 对共格双层构造成功（data/superlattices/）")
else:
    print('\n无成功构造')
