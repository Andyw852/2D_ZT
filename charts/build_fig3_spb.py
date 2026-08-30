# -*- coding: utf-8 -*-
"""SPB derivation, r=-1/2 acoustic-phonon: S = (kB/e)[2 F1(eta)/F0(eta) - eta]."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mplcache')
import math, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad

kB, e, h, L = 1.380649e-23, 1.602176634e-19, 6.62607015e-34, 2.44e-8
m0 = 9.109e-31

def F0(eta): return np.log1p(np.exp(np.asarray(eta, dtype=float)))
def F1(eta):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    out = np.zeros_like(eta)
    for i, et in enumerate(eta):
        v, _ = quad(lambda x: x/(np.exp(x-et)+1.0), 0, 80, limit=300)
        out[i] = v/math.gamma(2.0)
    return out if out.size > 1 else out[0]

def spb(eta, mstar, T, mu0=2.5e-3, mref=0.56, Tref=600.0):
    Nv = 2.0*(2*np.pi*mstar*m0*kB*T/h**2)**1.5           # m^-3
    n_m3 = Nv*F0(eta)
    S = -(kB/e)*(2.0*F1(eta)/F0(eta) - eta)               # electrons: negative
    mu = mu0*(mref/mstar)**2.5*(Tref/T)**1.5
    sg = n_m3*e*mu
    return n_m3/1e6, S*1e6, sg

def zt(eta, mstar, T, kL, mu0=2.5e-3):
    n, S, sg = spb(eta, mstar, T, mu0=mu0)
    ke = L*sg*T
    return n, S, sg, (S*1e-6)**2*sg*T/(ke+kL)

eta_scan = np.linspace(-5, 7, 160)
m_star = 0.8; kL0 = 1.5; T0 = 600.0
n600, S600, sg600, zt600 = zt(eta_scan, m_star, T0, kL0)
PF600 = (S600*1e-6)**2*sg600
opt_i = np.argmax(zt600); opt_n = n600[opt_i]

fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
axs = axs.ravel()
ax0 = axs[0]
ax0.plot(n600, np.abs(S600), color='#1f6fd6', lw=1.8, label='|S| (left)')
ax0.set_xscale('log'); ax0.set_xlabel('carrier concentration n (cm$^{-3}$)')
ax0.set_ylabel('|S|  (μV/K)', color='#1f6fd6'); ax0.tick_params(axis='y', labelcolor='#1f6fd6')
ax0b = ax0.twinx()
ax0b.plot(n600, zt600, color='#d64541', lw=1.8, label='ZT (right)')
ax0b.set_ylabel('ZT (κ$_L$=1.5)', color='#d64541'); ax0b.tick_params(axis='y', labelcolor='#d64541')
ax0.axvline(1e20, ls='--', lw=1, color='grey')
ax0.axvline(opt_n, ls=':', lw=1.2, color='k')
ax0.set_xlim(2e17, 3e21); ax0.set_ylim(0, 560)
ax0.text(1e20, 520, 'JARVIS n=1e20', fontsize=7.5, ha='center')
ax0.text(opt_n*1.15, 420, f'n*={opt_n:.1e}\ncm$^{{-3}}$', fontsize=8.5)
ax0.set_title(f'(a) |S| and ZT vs n (SPB Fermi integrals, acoustic phonon r=−½), T=600 K, m*={m_star} m$_e$')

ax1 = axs[1]
ax1.plot(n600, PF600, color='#2e8b57', lw=1.8, label='PF (left)')
ax1.set_xscale('log'); ax1.set_xlabel('carrier concentration n (cm$^{-3}$)')
ax1.set_ylabel('PF  (W m$^{-1}$ K$^{-2}$)', color='#2e8b57'); ax1.tick_params(axis='y', labelcolor='#2e8b57')
ax1b = ax1.twinx()
ax1b.plot(n600, sg600, color='#8a6d3b', lw=1.8, ls='--', label='σ (right)')
ax1b.set_yscale('log'); ax1b.set_ylabel('σ  (S/m)', color='#8a6d3b'); ax1b.tick_params(axis='y', labelcolor='#8a6d3b')
ax1.axvline(1e20, ls='--', lw=1, color='grey'); ax1.axvline(opt_n, ls=':', lw=1.2, color='k')
ax1.set_xlim(2e17, 3e21)
ax1.set_title('(b) PF = S²σ peaks at degenerate onset (~1e19–1e20 cm$^{-3}$)')

for T in [300, 450, 600, 750, 900]:
    nT, _, _, ztT = zt(eta_scan, m_star, T, kL0)
    axs[2].plot(nT, ztT, lw=1.6, label=f'T={T} K')
axs[2].set_xscale('log'); axs[2].set_xlabel('carrier concentration n (cm$^{-3}$)')
axs[2].set_ylabel('ZT'); axs[2].legend(frameon=False, fontsize=8)
axs[2].axvline(1e20, ls='--', lw=1, color='grey')
axs[2].set_xlim(2e17, 3e21)
axs[2].set_title('(c) optimum n* shifts up with T (higher T → higher optimal doping)')

T_scan = np.linspace(300, 1000, 90)
for kL in [0.3, 0.5, 1.0, 1.5, 3.0, 6.0]:
    ztT = np.array([zt(eta_scan, m_star, T, kL)[3].max() for T in T_scan])
    axs[3].plot(T_scan, ztT, lw=1.6, label=f'κ$_L$={kL}')
axs[3].set_xlabel('Temperature (K)'); axs[3].set_ylabel('max ZT over n')
axs[3].legend(frameon=False, fontsize=8, ncol=2)
axs[3].set_title('(d) ZT ceiling vs T: lowering κ$_L$ is the main design lever')

fig.suptitle('Numerical derivation: single-parabolic-band model (acoustic-phonon scattering)', fontsize=11.5, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig('charts/fig3_SPB_numerical_derivation.png', bbox_inches='tight')
plt.close(fig)

# print key numbers for the report
print(f'opt_n(600K) = {opt_n:.2e} cm-3, ZTmax = {zt600[opt_i]:.3f}')
for target in [1e18, 1e19, 3e19, 1e20, 3e20, 1e21]:
    j = np.argmin(np.abs(n600-target))
    print(f'n={target:.0e}: |S|={abs(S600[j]):.0f} uV/K, sigma={sg600[j]:.2e} S/m, PF={PF600[j]*1e3:.1f} mW/mK2, ZT={zt600[j]:.2f}')
