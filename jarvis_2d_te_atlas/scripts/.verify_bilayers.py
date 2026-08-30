#!/usr/bin/env python3
"""验证双层 POSCAR：原子间距、层间距、对称性。"""
import numpy as np, pandas as pd, glob
from pathlib import Path
from ase.io import read
from ase.neighborlist import neighbor_list
import spglib

OUT = Path('/home/wangchao/work_wc/2D_ZT/jarvis_2d_te_atlas/data/superlattices')
print(f'{"name":<28} {"nat":>3} {"min_d":>6} {"gap":>6} {"c":>7}  sym')
for f in sorted(glob.glob(str(OUT / '*_bilayer.POSCAR'))):
    a = read(f)
    # 最小原子间距（含层间）
    i, j, d = neighbor_list('ijd', a, cutoff=5.0)
    min_d = d.min() if len(d) else float('nan')
    # 跨层间隙
    z = a.positions[:, 2]
    zA = z[z < (z.min() + (z.max()-z.min())/2)]
    zB = z[z >= (z.min() + (z.max()-z.min())/2)]
    gap = zB.min() - zA.max() if len(zA) and len(zB) else float('nan')
    # 对称性（平面内周期，非周期 c）
    cell = a.cell[:2].tolist() + [[0,0,30.0]]
    cell = np.array(cell)
    pos = a.positions.copy(); pos[:, 2] %= 30.0
    sg = spglib.get_spacegroup((cell, pos, a.numbers), symprec=0.5)
    print(f'{f.split("/")[-1]:<28} {len(a):>3} {min_d:>6.2f} {gap:>6.2f} {a.cell[2][2]:>7.2f}  {sg}')
