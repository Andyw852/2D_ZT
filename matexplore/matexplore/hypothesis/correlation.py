# -*- coding: utf-8 -*-
"""假设核心：易算特征 <-> 热电目标 的相关面板 + 二维折线相关图。

这就是"把热电优质参数与材料本身的结构/电子/易算特征关联起来"的落地点：
  - 对每个 (特征, 目标, carrier) 计算 Spearman 相关
  - 输出 CSV + 二维折线相关图(横轴=特征，纵轴=Spearman rho，多折线=目标/载流子)
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# 特征面板(与 scripts/build_feature_target_analysis.py 一致)
FEATURES = [
    ("m* (median, m_e)", "m_star"),
    ("Eg (band gap, eV)", "Eg_optb88vdw"),
    ("electronegativity mean", "electronegativity_mean"),
    ("electronegativity range", "electronegativity_range"),
    ("atomic mass mean", "atomic_mass_mean"),
    ("atomic mass max", "atomic_mass_max"),
    ("atomic radius mean", "atomic_radius_mean"),
    ("Z (atomic number) mean", "Z_mean"),
    ("ionization energy mean", "ionization_energy_mean"),
    ("electron affinity mean", "electron_affinity_mean"),
    ("row mean", "row_mean"),
    ("group range", "group_range"),
    ("density (g/cm3)", "density"),
    ("formation energy /atom (eV)", "formation_peratom"),
    ("exfoliation energy (eV)", "exfoliation"),
    ("n_sites", "nsites"),
    ("space group number", "spg_number"),
    ("in-plane area (A^2)", "inplane_area"),
    ("vacuum / sqrt(area)", "aspect_c_over_sqrtA"),
]

TARGETS = [
    ("ZT_e", "ZT_e_{c}"),
    ("log10 PF", "PF_mean_{c}"),
    ("|S| (uV/K)", "S_median_{c}"),
    ("log10 sigma", "log_sigma_dom_geo_{c}"),
]


def compute_correlation(panel):
    rows = []
    for carrier in ["n", "p"]:
        d = panel.copy()
        d["m_star"] = d["m_elec_median"] if carrier == "n" else d["m_hole_median"]
        for feat_label, feat_col in FEATURES:
            for tlabel, tcol in TARGETS:
                tcol = tcol.format(c=carrier)
                sub = d[[feat_col, tcol]].dropna()
                if len(sub) < 50:
                    continue
                x = sub[feat_col].to_numpy(float)
                y = sub[tcol].to_numpy(float)
                if tlabel == "log10 PF":
                    y = np.log10(np.clip(y, 1e-9, None))
                elif tlabel == "|S| (uV/K)":
                    y = np.abs(y)
                rho = spearmanr(x, y).statistic
                if np.isfinite(rho):
                    rows.append(dict(carrier=carrier, feature=feat_label,
                                     target=tlabel, spearman=float(rho), n=len(sub)))
    return pd.DataFrame(rows)


def plot_polyline(corr, out_png):
    feat_order = (corr.groupby("feature")["spearman"]
                  .apply(lambda s: s.abs().max()).sort_values(ascending=False).index.tolist())
    carrier_color = {"n": "#1f77b4", "p": "#d62728"}
    carrier_ls = {"n": "-", "p": "--"}
    targets = ["ZT_e", "log10 PF", "|S| (uV/K)", "log10 sigma"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    for ax, tgt in zip(axes.ravel(), targets):
        sub = corr[corr.target == tgt]
        for c in ["n", "p"]:
            cc = sub[sub.carrier == c].set_index("feature")["spearman"].reindex(feat_order)
            ax.plot(range(len(feat_order)), cc.to_numpy(), ls=carrier_ls[c],
                    color=carrier_color[c], marker="o", ms=4, lw=1.4, label=f"{c}-type")
        ax.axhline(0, color="gray", lw=0.8)
        ax.axhline(0.5, color="gray", lw=0.6, ls=":")
        ax.axhline(-0.5, color="gray", lw=0.6, ls=":")
        ax.set_xticks(range(len(feat_order)))
        ax.set_xticklabels([f.split("(")[0].strip() for f in feat_order],
                           rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Spearman rho")
        ax.set_title(f"target: {tgt}")
        ax.set_ylim(-1.05, 1.05)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)
    fig.suptitle("Easy-feature vs thermoelectric-target Spearman correlation "
                 "(JARVIS dft_2d, 1103 2D materials)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def form_hypothesis(corr, top_n=8):
    """从相关面板提炼"假设"：每个目标给出最相关的若干特征及其方向。

    返回 dict: { target: [ {feature, carrier, rho}, ... ] }
    """
    hyp = {}
    for tgt in TARGETS:
        tlabel = tgt[0]
        sub = corr[corr.target == tlabel].copy()
        sub["absrho"] = sub.spearman.abs()
        top = sub.sort_values("absrho", ascending=False).head(top_n)
        hyp[tlabel] = top[["feature", "carrier", "spearman", "n"]].to_dict("records")
    return hyp
