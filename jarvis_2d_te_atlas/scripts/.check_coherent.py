#!/usr/bin/env python3
import pandas as pd
for mod in ('ase', 'pymatgen', 'spglib'):
    try:
        m = __import__(mod)
        print(mod, 'OK', getattr(m, '__version__', ''))
    except Exception as e:
        print(mod, 'MISSING', e)
df = pd.read_csv('data/processed/superlattice_candidate_pairs.csv')
print('总对:', len(df))
c = df[(df.lattice_mismatch_pct < 1) & (df.angle_mismatch_deg < 1)]
print('晶格失配<1% 且 角度失配<1度 的共格对:', len(c))
cols = ['A_jid','A_formula','B_jid','B_formula','lattice_mismatch_pct','area_mismatch_pct','angle_mismatch_deg','S_contrast','A_ZT_e','B_ZT_e']
print(c[cols].head(15).to_string(index=False))
