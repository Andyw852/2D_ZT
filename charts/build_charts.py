# -*- coding: utf-8 -*-
"""Build the full chart suite answering the 5 design questions.
All figures use English labels (no CJK font available); Chinese explanation lives in the report.
"""
import os, sys, ast, warnings
warnings.filterwarnings('ignore')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mplcache')
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
from scipy.stats import spearmanr

ROOT = 'jarvis_2d_te_atlas'
OUT = 'charts'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9.5,
                     'figure.dpi': 130, 'savefig.dpi': 130, 'axes.spines.top': False,
                     'axes.spines.right': False})

z = pd.read_csv(f'{ROOT}/data/processed/ZT_e_all.csv')
z['m_use'] = z.apply(lambda r: r.m_elec_median if r.carrier == 'n' else r.m_hole_median, axis=1)
zn = z[z.carrier == 'n'].copy(); zp = z[z.carrier == 'p'].copy()

# ---------------- Fig 1: necessary electronic conditions (Q1) ----------------
fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.5))
ax = ax.ravel()
# (a) ZT_e vs Eg
for df, c, lab in [(zn, '#1f6fd6', 'n-type'), (zp, '#d64541', 'p-type')]:
    ax[0].scatter(df.Eg_optb88vdw, df.ZT_e, s=7, alpha=0.35, c=c, label=lab, edgecolors='none')
ax[0].axvspan(1.0, 2.0, color='gold', alpha=0.18)
ax[0].axvline(1.0, ls='--', lw=0.8, color='k', alpha=0.5); ax[0].axvline(2.0, ls='--', lw=0.8, color='k', alpha=0.5)
ax[0].set_xlabel('Band gap  E$_g$  (eV, OptB88vdW)'); ax[0].set_ylabel('ZT$_e$  = S$^2$/L  (electronic ceiling)')
ax[0].set_ylim(0, 30)
ax[0].text(1.5, 28.0, 'sweet spot\nEg ≈ 1–2 eV', ha='center', fontsize=8.5, color='#7a5c00')
ax[0].legend(frameon=False, loc='upper right', fontsize=8)
ax[0].set_title('(a) High ZT$_e$ requires semiconducting gap  (Spearman n +0.63 / p +0.76)')
# (b) ZT_e vs m*
for df, c, lab in [(zn, '#1f6fd6', 'n-type'), (zp, '#d64541', 'p-type')]:
    dd = df[df.m_use > 0]
    ax[1].scatter(dd.m_use, dd.ZT_e, s=7, alpha=0.35, c=c, label=lab, edgecolors='none')
ax[1].set_xscale('log'); ax[1].set_xlim(0.05, 60)
ax[1].set_xlabel('Effective mass  m*  (m$_e$, median of principal values)')
ax[1].set_ylabel('ZT$_e$')
ax[1].set_ylim(0, 30); ax[1].legend(frameon=False, loc='upper left', fontsize=8)
ax[1].set_title('(b) Heavier m* raises ZT$_e$ but kills σ  (Spearman +0.78 / σ −0.87~−0.90)')
# (c) S-sigma tradeoff
zz = pd.concat([zn, zp])
sc = ax[2].scatter(zz['log_sigma_dom_geo'], zz.S_median.abs(), c=np.log10(zz.ZT_e+1e-6), s=8,
                   cmap='magma', alpha=0.75, edgecolors='none')
cb = plt.colorbar(sc, ax=ax[2]); cb.set_label('log$_{10}$ ZT$_e$')
ax[2].set_xlabel('log$_{10}$ σ  (S/m, dominant-channel)')
ax[2].set_ylabel('|S|  (μV/K)')
ax[2].set_yscale('log')
ax[2].set_title('(c) S–σ trade-off: high ZT$_e$ = high |S| + medium σ (diagonal band)')
# (d) PF vs Eg
for df, c, lab in [(zn, '#1f6fd6', 'n-type'), (zp, '#d64541', 'p-type')]:
    ax[3].scatter(df.Eg_optb88vdw, df.PF_mean, s=7, alpha=0.35, c=c, label=lab, edgecolors='none')
ax[3].axvspan(1.0, 2.0, color='gold', alpha=0.18)
ax[3].set_yscale('log'); ax[3].set_xlabel('Band gap  E$_g$  (eV)')
ax[3].set_ylabel('PF  (μW m$^{-1}$ K$^{-2}$)')
ax[3].legend(frameon=False, loc='upper right', fontsize=8)
ax[3].set_title('(d) PF peaks at intermediate gap (S·σ product, Spearman +0.45/+0.49)')
fig.tight_layout()
fig.savefig(f'{OUT}/fig1_ZT_conditions_sweetspot.png', bbox_inches='tight')
plt.close(fig)
print('fig1 done')

# ---------------- Fig 2: feature→target correlation polylines (Q2) ----------------
fc = pd.read_csv('matexplore/reports/feature_target_correlation.csv')
cat = {'Eg (band gap, eV)': 'electronic', 'm* (median, m_e)': 'electronic',
       'electronegativity mean': 'composition', 'electronegativity range': 'composition',
       'atomic mass mean': 'composition', 'atomic mass max': 'composition',
       'atomic radius mean': 'composition', 'Z (atomic number) mean': 'composition',
       'ionization energy mean': 'composition', 'electron affinity mean': 'composition',
       'row mean': 'composition', 'group range': 'composition',
       'density (g/cm3)': 'composition', 'formation energy /atom (eV)': 'composition',
       'exfoliation energy (eV)': 'geometry', 'n_sites': 'geometry',
       'space group number': 'geometry', 'in-plane area (Å^2)': 'geometry',
       'vacuum / sqrt(area)': 'geometry'}
fc['cat'] = fc.feature.map(cat)
order = fc[(fc.target == 'ZT_e') & (fc.carrier == 'n')].sort_values('spearman').feature.tolist()
targets = ['ZT_e', '|S| (uV/K)', 'log10 PF', 'log10 sigma']
titles = ['ZT$_e$  (electronic ceiling)', '|S|  (μV/K)', 'log$_{10}$ PF', 'log$_{10}$ σ']
fig, axs = plt.subplots(2, 2, figsize=(13, 9))
axs = axs.ravel()
band = {'electronic': ('#c9e4ff', 'electronic: Eg, m*'),
        'composition': ('#ffe9c9', 'composition (Magpie, no DFT)'),
        'geometry': ('#d5f5d5', 'geometry / exfoliation')}
for ax, tgt, tt in zip(axs, targets, titles):
    d = fc[fc.target == tgt].pivot_table(index='feature', columns='carrier', values='spearman').reindex(order)
    x = np.arange(len(order))
    for i, f in enumerate(order):
        c, lab = band[cat[f]]
        ax.axvspan(i-0.5, i+0.5, color=c, alpha=0.45)
    ax.axhline(0, color='k', lw=0.8)
    ax.plot(x, d['n'], '-o', ms=3.5, lw=1.6, color='#1f6fd6', label='n-type')
    ax.plot(x, d['p'], '--s', ms=3.5, lw=1.6, color='#d64541', label='p-type')
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=55, ha='right', fontsize=7.5)
    ax.set_ylim(-1, 1); ax.set_ylabel('Spearman ρ')
    ax.set_title(tt, fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc='lower left')
handles = [plt.Line2D([0], [0], color=c, lw=6, alpha=0.5, label=lab) for c, lab in band.values()]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
fig.suptitle('Feature → thermoelectric target correlations (Spearman, 1103 2D materials; blue=n, red=p)', y=1.0, fontsize=11)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
fig.savefig(f'{OUT}/fig2_feature_target_polyline.png', bbox_inches='tight')
plt.close(fig)
print('fig2 done')

# ---------------- Fig 4: anisotropy + lattice-preferred pairing (Q5) ----------------
an = pd.read_csv(f'{ROOT}/data/audit/conductivity_spectrum_audit_n.csv', usecols=['jid','ncond_D','ncond_A_dom'])
ap = pd.read_csv(f'{ROOT}/data/audit/conductivity_spectrum_audit_p.csv', usecols=['jid','pcond_D','pcond_A_dom'])
zn2 = zn.merge(an, on='jid', how='inner').dropna(subset=['ncond_D','ncond_A_dom'])
zp2 = zp.merge(ap, on='jid', how='inner').dropna(subset=['pcond_D','pcond_A_dom'])
rho_n = spearmanr(zn2.ZT_e, zn2.ncond_D)[0]; rho_p = spearmanr(zp2.ZT_e, zp2.pcond_D)[0]
sl = pd.read_csv(f'{ROOT}/data/processed/superlattice_candidate_pairs.csv')
slc = sl[(sl.lattice_mismatch_pct < 1) & (sl.area_mismatch_pct < 1)]
fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
axs = axs.ravel()
# (a) D vs A_dom
axs[0].scatter(zn2.ncond_D, zn2.ncond_A_dom, s=8, alpha=0.5, c='#1f6fd6', label='n-type', edgecolors='none')
axs[0].scatter(zp2.pcond_D, zp2.pcond_A_dom, s=8, alpha=0.5, c='#d64541', label='p-type', edgecolors='none')
axs[0].set_xlabel('suppressed-channel contrast  D = log$_{10}$(σ$_2$/σ$_1$)')
axs[0].set_ylabel('dominant-channel anisotropy  A$_{dom}$ = log$_{10}$(σ$_3$/σ$_2$)')
axs[0].set_xlim(-0.2, 6.5); axs[0].axhline(0.3, ls='--', lw=0.8, color='grey')
axs[0].text(4.6, 2.6, 'quasi-2D corner:\none channel suppressed,\ntwo dominant ≈ isotropic', fontsize=8)
axs[0].legend(frameon=False, fontsize=8, loc='upper left')
axs[0].set_title('(a) Conductivity principal spectrum: quasi-2D signature (D≈3.1–3.3, A$_{dom}$≈0)')
# (b) ZT_e vs D
axs[1].scatter(zn2.ncond_D, zn2.ZT_e, s=8, alpha=0.45, c='#1f6fd6', label='n-type', edgecolors='none')
axs[1].scatter(zp2.pcond_D, zp2.ZT_e, s=8, alpha=0.45, c='#d64541', label='p-type', edgecolors='none')
axs[1].set_xlabel('suppressed-channel contrast D')
axs[1].set_ylabel('ZT$_e$'); axs[1].set_ylim(0, 30)
axs[1].text(0.05, 0.95, f'Spearman\nn {rho_n:+.2f} / p {rho_p:+.2f}', transform=axs[1].transAxes, fontsize=9,
            va='top', bbox=dict(fc='white', ec='grey', alpha=0.8))
axs[1].legend(frameon=False, fontsize=8, loc='upper right')
axs[1].set_title('(b) Is strong anisotropy needed for high ZT$_e$? (weak)')
# (c) lattice-parameter prefiltered pairs（晶格参数预筛，非已构造超晶格）
axs[2].scatter(sl.lattice_mismatch_pct, sl.S_contrast, s=8, alpha=0.5, c='#888888', label=f'all {len(sl)} pairs', edgecolors='none')
axs[2].scatter(slc.lattice_mismatch_pct, slc.S_contrast, s=14, alpha=0.85, c='#2e8b57', label='coherent (lattice+area <1%)', edgecolors='k', linewidths=0.3)
axs[2].axvline(1.0, ls='--', lw=1, color='#2e8b57')
axs[2].set_xlabel('lattice mismatch (%)')
axs[2].set_ylabel('|S$_n$ − S$_p$| contrast  (μV/K)')
axs[2].legend(frameon=False, fontsize=8, loc='lower left')
axs[2].set_title(f'(c) Lattice-preferred pairing pool: {len(slc)} coherent pairs with big S contrast')
# (d) bipolar materials
piv = z.pivot_table(index='jid', columns='carrier', values='ZT_e').dropna(subset=['n','p'])
axs[3].scatter(piv['n'], piv['p'], s=8, alpha=0.45, c='#6b6bd6', edgecolors='none')
axs[3].axline((0, 0), slope=1, ls='--', lw=0.8, color='grey')
bip = piv[(piv.n > 2) & (piv.p > 2)]
axs[3].scatter(bip.n, bip.p, s=16, alpha=0.9, c='#2e8b57', edgecolors='k', linewidths=0.3, label=f'n & p both ZT$_e$>2 ({len(bip)})')
axs[3].set_xlabel('n-type ZT$_e$'); axs[3].set_ylabel('p-type ZT$_e$'); axs[3].set_xlim(0, 30); axs[3].set_ylim(0, 30)
axs[3].legend(frameon=False, fontsize=8, loc='upper left')
axs[3].set_title('(d) Bipolar materials: one material, both carrier types — tunability seed')
fig.suptitle('Anisotropy landscape & lattice-preferred pairing pool', fontsize=11.5, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f'{OUT}/fig4_anisotropy_superlattice.png', bbox_inches='tight')
plt.close(fig)
print('fig4 done; coherent pairs =', len(slc))

# ---------------- Fig 5: kappa_L design map (Q3) ----------------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={'width_ratios': [1.15, 1]})
ax = axs[0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.text(5, 9.6, 'κ$_L$ = ⅓ Σ C$_v$·v$_g$²·τ   →   what a 2D design can control', ha='center', fontsize=11, weight='bold')
boxes = [
    (0.2, 6.6, 'SHAPE\nwrinkle / planar / skeleton\n→ v$_g$ (dispersion), 2D confinement', '#c9e4ff'),
    (3.5, 6.6, 'HOMOGENEITY\nporous? soft building blocks?\n→ τ (scattering), rattlers', '#ffe9c9'),
    (6.8, 6.6, 'BONDING\nvibration frequency, acoustic /\noptical branches, ω$_D$\n→ C$_v$, v$_g$', '#d5f5d5')]
for x, y, t, c in boxes:
    ax.add_patch(FancyBboxPatch((x, y), 3.0, 2.1, boxstyle='round,pad=0.05', fc=c, ec='grey', lw=0.8))
    ax.text(x+1.5, y+1.05, t, ha='center', va='center', fontsize=8.2)
ax.annotate('', xy=(1.7, 6.5), xytext=(1.7, 5.9), arrowprops=dict(arrowstyle='->', lw=1.4))
ax.annotate('', xy=(5.0, 6.5), xytext=(5.0, 5.9), arrowprops=dict(arrowstyle='->', lw=1.4))
ax.annotate('', xy=(8.3, 6.5), xytext=(8.3, 5.9), arrowprops=dict(arrowstyle='->', lw=1.4))
ax.add_patch(FancyBboxPatch((1.2, 3.4), 7.6, 2.2, boxstyle='round,pad=0.05', fc='#f5f5f5', ec='grey', lw=0.8))
ax.text(5, 4.5, 'lower κ$_L$  →  higher ZT at fixed PF\n(low v$_g$: heavy atoms, soft modes;  short τ: porosity, interfaces, anharmonicity;\nsmall C$_v$: low ω$_D$ — but optical branches add channels: care needed)', ha='center', va='center', fontsize=8.2)
ax.text(1.2, 2.6, '⚠ not in JARVIS dft_2d (κ$_L$=0 coverage). Sources: phonon calc (phonopy/DFPT) or experimental κ$_L$ curves (Starrydata2, 6740 curves).', fontsize=8, color='#8a1f1f')
ax.set_title('(a) Phonon-side design map (physics formula, data-gap marked)', fontsize=10.5)
# (b) experimental kappa_L curves
import ast as _ast
cur = pd.read_csv(f'{ROOT}/data/raw/external/starrydata2/ThermoelectricMaterials_curves.csv',
                  usecols=['composition','prop_x','prop_y','x','y'])
kcur = cur[(cur.prop_x == 'Temperature') & (cur.prop_y == 'Lattice thermal conductivity')]
keys = [('Bi2Te3', '#d64541'), ('PbTe', '#1f6fd6'), ('SnSe', '#2e8b57'), ('CoSb3', '#8a6d3b')]
plotted = {}
for comp, c in keys:
    m = kcur[kcur.composition.str.contains(comp, case=False, na=False, regex=False)]
    cnt = 0
    for _, row in m.iterrows():
        try:
            x = _ast.literal_eval(row.x); y = _ast.literal_eval(row.y)
        except Exception:
            continue
        x = np.array(x, dtype=float); y = np.array(y, dtype=float)
        ok = (x > 250) & (x < 1000) & (y > 0) & (y < 30)
        if ok.sum() < 5: continue
        axs[1].plot(x[ok], y[ok], lw=1.0, alpha=0.6, c=c)
        cnt += 1
        if cnt >= 3: break
    plotted[comp] = cnt
axs[1].set_xlabel('T (K)'); axs[1].set_ylabel('κ$_L$ (W m$^{-1}$ K$^{-1}$)')
axs[1].set_ylim(0, 15)
for comp, c in keys:
    axs[1].plot([], [], lw=1.6, c=c, label=f'{comp} ({plotted[comp]} curves)')
axs[1].legend(frameon=False, fontsize=8)
axs[1].set_title('(b) Experimental κ$_L$(T) from Starrydata2 — the missing anchor')
fig.tight_layout()
fig.savefig(f'{OUT}/fig5_kL_design_map.png', bbox_inches='tight')
plt.close(fig)
print('fig5 done')
print('ALL FIGURES WRITTEN to', OUT)
