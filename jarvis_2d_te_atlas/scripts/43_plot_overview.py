"""美观版设计规律图（2×2）+ 独立高 ZT_e 流形主图。"""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from pathlib import Path

root = Path(__file__).resolve().parents[1]
figdir = root / 'figures'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.linewidth': 0.8, 'axes.edgecolor': '#bbbbbb',
    'xtick.color': '#444444', 'ytick.color': '#444444',
    'axes.labelsize': 12, 'axes.titlesize': 12.5, 'legend.frameon': False,
})

n = pd.read_parquet(root/'manifolds/n_atlas_consensus_zte.parquet')
p = pd.read_parquet(root/'manifolds/p_atlas_consensus_zte.parquet')
sp = pd.read_csv(root/'data/processed/superlattice_candidate_pairs.csv')
zt = pd.concat([n.assign(carrier='n'), p.assign(carrier='p')], ignore_index=True)
zt['logZT'] = np.log10(np.clip(zt['ZT_e'], 1e-3, None))

CMAP = 'magma'
NORM = mcolors.Normalize(vmin=-2.5, vmax=1.2)

def despine(ax):
    for s in ['top','right']:
        ax.spines[s].set_visible(False)
    ax.grid(True, ls=':', lw=0.5, color='#e0e0e0', alpha=0.7)
    ax.set_axisbelow(True)

# ============ 图 1: 2×2 总览 ============
fig, axes = plt.subplots(2, 2, figsize=(13, 11.5))

# (a) p-type manifold
ax = axes[0,0]; despine(ax)
sc = ax.scatter(p['Phi_1'], p['Phi_2'], c=np.log10(np.clip(p['ZT_e'],1e-3,None)),
                cmap=CMAP, norm=NORM, s=16, alpha=0.82, linewidths=0, rasterized=True)
hi = p[p['ZT_e']>2]
if len(hi):
    cx, cy = hi['Phi_1'].mean(), hi['Phi_2'].mean()
    ax.scatter([cx],[cy], marker='*', s=300, color='#ffd166', edgecolor='#333', linewidth=1, zorder=6)
    ax.annotate('high-ZT$_{e}$\nregion', xy=(cx,cy), xytext=(cx+0.045, cy+0.02),
                fontsize=10.5, fontweight='bold', color='#7a4a00',
                arrowprops=dict(arrowstyle='->', color='#7a4a00', lw=1.3))
ax.set_xlabel(r'$\Phi_1$  (electronic-structure axis)'); ax.set_ylabel(r'$\Phi_2$')
ax.set_title('(a)  p-type joint manifold, coloured by ZT$_{e}$')
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r'log$_{10}$ ZT$_{e}$', fontsize=10)
cb.ax.tick_params(labelsize=9)

# (b) ZT_e vs Eg
ax = axes[0,1]; despine(ax)
for c, col, lab in [('n','#3b6ea5','n-type'), ('p','#c44536','p-type')]:
    sub = zt[zt['carrier']==c]
    ax.scatter(sub['Eg_optb88vdw'], sub['ZT_e'], c=col, s=14, alpha=0.6, linewidths=0, label=lab, rasterized=True)
ax.set_yscale('log'); ax.set_ylim(1e-3, 40)
ax.set_xlabel(r'band gap  $E_g$ (eV)'); ax.set_ylabel(r'ZT$_{e}$ ceiling (log)')
ax.set_title('(b)  high ZT$_{e}$ requires a band gap')
ax.axvline(0.05, color='#999999', ls=':', lw=1)
ax.text(0.12, 20, 'metals\n$(S\\approx 0)$', fontsize=9, color='#888888')
ax.legend(loc='upper left', fontsize=10, markerscale=1.6)

# (c) S-σ trade-off
ax = axes[1,0]; despine(ax)
sc2 = ax.scatter(zt['absS'], zt['log_sigma_dom_geo'], c=zt['logZT'],
                 cmap=CMAP, norm=NORM, s=16, alpha=0.78, linewidths=0, rasterized=True)
ax.set_xlabel(r'$|S|$  (μV/K)'); ax.set_ylabel(r'log$_{10}$ conductivity scale')
ax.set_title('(c)  Seebeck–conductivity trade-off')
cb2 = fig.colorbar(sc2, ax=ax, pad=0.02); cb2.set_label(r'log$_{10}$ ZT$_{e}$', fontsize=10)
cb2.ax.tick_params(labelsize=9)

# (d) superlattice candidates
ax = axes[1,1]; despine(ax)
sp['pair_ZT'] = np.minimum(sp['A_ZT_e'], sp['B_ZT_e'])
sc3 = ax.scatter(sp['lattice_mismatch_pct'], sp['S_contrast'], c=sp['pair_ZT'],
                 cmap='viridis', s=34, alpha=0.9, linewidths=0.6, edgecolor='white', rasterized=True)
ax.set_xlabel('lattice mismatch (%)'); ax.set_ylabel(r'$|S_n - S_p|$ (μV/K)')
ax.set_title('(d)  lattice-parameter prefiltered pairs  (n = 89)')
ax.axvline(0.5, color='#c44536', ls='--', lw=1, alpha=0.8)
ax.text(0.55, 545, 'mismatch < 0.5%', fontsize=9, color='#c44536')
cb3 = fig.colorbar(sc3, ax=ax, pad=0.02); cb3.set_label('min(ZT$_{e}$)', fontsize=10)
cb3.ax.tick_params(labelsize=9)

fig.suptitle('Design rules for 2D thermoelectrics — JARVIS dft_2d (N = 1103)',
             fontsize=15, fontweight='bold', y=0.995)
fig.tight_layout(rect=[0,0,1,0.97])
for ext in ['png','pdf']:
    fig.savefig(figdir/f'design_rules_overview.{ext}', dpi=220, bbox_inches='tight')
plt.close(fig)
print('saved design_rules_overview')

# ============ 图 2: 独立高 ZT_e 流形主图（n/p 并排）============
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharey=False)
for ax, df, lab, letter in [(axes[0], n, 'n-type', '(a)'), (axes[1], p, 'p-type', '(b)')]:
    despine(ax)
    sc = ax.scatter(df['Phi_1'], df['Phi_2'], c=np.log10(np.clip(df['ZT_e'],1e-3,None)),
                    cmap=CMAP, norm=NORM, s=20, alpha=0.85, linewidths=0, rasterized=True)
    hi = df[df['ZT_e']>2]
    if len(hi):
        cx, cy = hi['Phi_1'].mean(), hi['Phi_2'].mean()
        ax.scatter([cx],[cy], marker='*', s=340, color='#ffd166', edgecolor='#333', linewidth=1.2, zorder=6)
    ax.set_xlabel(r'$\Phi_1$'); ax.set_ylabel(r'$\Phi_2$')
    ax.set_title(f'{letter}  {lab} joint manifold  ·  coloured by ZT$_{{e}}$')
fig.suptitle('High-ZT$_{e}$ manifold: the electronic-structure axis separates the thermoelectric corner',
             fontsize=13.5, fontweight='bold', y=1.02)
cb = fig.colorbar(sc, ax=axes, pad=0.02, aspect=40); cb.set_label(r'log$_{10}$ ZT$_{e}$ ceiling', fontsize=11)
fig.tight_layout(rect=[0,0,1,0.99])
for ext in ['png','pdf']:
    fig.savefig(figdir/f'high_ZT_manifold.{ext}', dpi=220, bbox_inches='tight')
plt.close(fig)
print('saved high_ZT_manifold')
