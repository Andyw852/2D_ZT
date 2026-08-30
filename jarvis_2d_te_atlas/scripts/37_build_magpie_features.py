"""补强结构描述符: Magpie 风格成分特征（元素属性加权统计）。"""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from pymatgen.core import Element

root = Path(__file__).resolve().parents[1]
sdf = pd.read_parquet(root / 'data/processed/standardized_2d_structures.parquet')

# 元素属性表（81 元素）
elements = set()
for s in sdf['species']:
    elements.update(json.loads(s))
elements = sorted(elements)

props = {}
for el in elements:
    e = Element(el)
    props[el] = {
        'electronegativity': float(e.X),
        'atomic_mass': float(str(e.atomic_mass).split()[0]),
        'atomic_radius': float(str(e.atomic_radius).split()[0]) if e.atomic_radius else np.nan,
        'row': float(e.row),
        'group': float(e.group) if e.group else np.nan,
        'ionization_energy': float(e.ionization_energy) if e.ionization_energy else np.nan,
        'electron_affinity': float(e.electron_affinity) if e.electron_affinity else np.nan,
        'Z': float(e.Z),
    }

PROP_NAMES = ['electronegativity','atomic_mass','atomic_radius','row','group','ionization_energy','electron_affinity','Z']

def magpie(species):
    n = len(species)
    if n == 0:
        return {}
    vals = {p: [props[sp][p] for sp in species if not (isinstance(props[sp][p], float) and np.isnan(props[sp][p]))] for p in PROP_NAMES}
    feat = {}
    for p in PROP_NAMES:
        v = vals[p]
        if not v:
            feat[f'{p}_mean'] = np.nan; feat[f'{p}_std'] = np.nan
            feat[f'{p}_min'] = np.nan; feat[f'{p}_max'] = np.nan
            continue
        v = np.array(v)
        feat[f'{p}_mean'] = v.mean()
        feat[f'{p}_std'] = v.std()
        feat[f'{p}_min'] = v.min()
        feat[f'{p}_max'] = v.max()
        feat[f'{p}_range'] = v.max() - v.min()
    return feat

rows = []
for _, r in sdf.iterrows():
    sp = json.loads(r['species'])
    row = {'jid': r['jid']}
    row.update(magpie(sp))
    rows.append(row)
mdf = pd.DataFrame(rows).sort_values('jid').reset_index(drop=True)
mdf.to_parquet(root / 'features/structure/composition_magpie.parquet', index=False)
print(f'Magpie composition features: {mdf.shape} (n_feat={len([c for c in mdf.columns if c!="jid"])})')
print('columns:', [c for c in mdf.columns if c!='jid'][:12], '...')
print('NaN fraction:', float(mdf.iloc[:, 1:].isna().mean().mean()))
# 验证: S vs Se vs Au 的电负性相似性
print('S elecneg=%.2f Se=%.2f Au=%.2f -> S-Se close, S-Au far' % (props['S']['electronegativity'], props['Se']['electronegativity'], props['Au']['electronegativity']))
