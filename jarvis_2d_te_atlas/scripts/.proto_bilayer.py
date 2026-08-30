#!/usr/bin/env python3
"""原型：从 snapshot 取一对材料，旋转+应变对齐晶格，构造共格双层，写 POSCAR。"""
import json, numpy as np, sys
from pathlib import Path
from ase import Atoms
from ase.io import write

ROOT = Path('/home/wangchao/work_wc/2D_ZT/jarvis_2d_te_atlas')
snap = json.load(open(ROOT/'data/raw/jarvis/dft_2d_snapshot.json'))
by_id = {e['id']: e['attributes'] for e in snap}

def get_structure(jid):
    att = by_id.get(jid)
    if att is None: return None
    lat = np.array(att['lattice_vectors'], dtype=float)
    cart = np.array(att['cartesian_site_positions'], dtype=float)
    species = att['species_at_sites']
    return lat, cart, species

def plane_aligned(lat):
    """返回 (a2,b2) 平面内基矢、真空轴 z, 以及是否平面内 z~0"""
    a, b, c = lat[0], lat[1], lat[2]
    plane_ok = abs(a[2]) < 1e-3 and abs(b[2]) < 1e-3
    z_ok = abs(c[0]) < 1e-3 and abs(c[1]) < 1e-3
    return plane_ok, z_ok

jidA, jidB = 'dft_2d_JVASP-5888', 'dft_2d_JVASP-77690'  # OSn / BrOSb
la, ca, sa = get_structure(jidA)
lb, cb, sb = get_structure(jidB)
print('A', jidA, sa, 'lat a,b,c len:', np.linalg.norm(la[0]), np.linalg.norm(la[1]), np.linalg.norm(la[2]))
print('B', jidB, sb, 'lat a,b,c len:', np.linalg.norm(lb[0]), np.linalg.norm(lb[1]), np.linalg.norm(lb[2]))
pa, za = plane_aligned(la); pb, zb = plane_aligned(lb)
print('A plane_aligned:', pa, za, '| B plane_aligned:', pb, zb)

# 平面内基矢 (2D)
a1, b1 = la[0][:2], la[1][:2]
a2, b2 = lb[0][:2], lb[1][:2]
ang = lambda v: np.degrees(np.arctan2(v[1], v[0]))
print('a1 角:', round(ang(a1),2), 'b1 角:', round(ang(b1),2), '| a2 角:', round(ang(a2),2), 'b2 角:', round(ang(b2),2))
print('|a1|', round(np.linalg.norm(a1),4), '|b1|', round(np.linalg.norm(b1),4), '| |a2|', round(np.linalg.norm(a2),4), '|b2|', round(np.linalg.norm(b2),4))

# 应变 S: L1 = S @ L2  (2x2), L2 -> L1
L1 = np.array([a1, b1]); L2 = np.array([a2, b2])
S = L1 @ np.linalg.inv(L2)
print('应变矩阵 S:'); print(np.round(S, 5))
print('|S-I| fro:', round(np.linalg.norm(S - np.eye(2)), 4))

# B 原子平面内应变 + 层叠
cb2 = cb.copy()
cb2[:, :2] = cb2[:, :2] @ S.T   # cart = S @ cart
# 层间距：B 的原子 z 抬升到 A 顶面之上 3.3 A
zA_top = ca[:, 2].max()
zB_bot = cb2[:, 2].min()
gap = 3.3
shift = zA_top + gap - zB_bot
cb2[:, 2] += shift
zB_top = cb2[:, 2].max()
# 真空层 15 A
c_len = zB_top + 15.0
new_cell = la.copy()
new_cell[2] = [0, 0, c_len]
all_pos = np.vstack([ca, cb2])
all_sp = list(sa) + list(sb)
atoms = Atoms(all_sp, positions=all_pos, cell=new_cell, pbc=[True, True, False])
print('双层原子数:', len(atoms), '晶格 c:', round(c_len,2))
print('层内最小跨层间距:', round(min(np.linalg.norm(p1-p2) for p1 in cb2 for p2 in ca), 3))
from ase.io import write
write('/home/wangchao/work_wc/2D_ZT/jarvis_2d_te_atlas/data/superlattices/OSn_BrOSb_bilayer.POSCAR', atoms, vasp5=True, direct=True)
print('POSCAR 已写')
