"""Run numerical scans, response surfaces and report generation."""

from __future__ import annotations

import json
import os
from dataclasses import asdict, replace
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/zt_deep_physics_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from zt_model import (
    ElectronicParams,
    LatticeParams,
    StructureParams,
    bipolar_gap_threshold_eV,
    electronic_transport,
    lattice_transport,
    optimise_eta,
    spb_seebeck_lorenz,
)


PROJECT = Path(__file__).resolve().parents[1]
REPO = PROJECT.parent
CONFIG = PROJECT / "config" / "baseline.json"
OUT = PROJECT / "outputs"
FIG = PROJECT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 180,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def load_params():
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    return (
        cfg,
        ElectronicParams(**cfg["electronic"]),
        LatticeParams(**cfg["lattice"]),
        StructureParams(**cfg["structure"]),
    )


def eta_grid(cfg) -> np.ndarray:
    scan = cfg["scan"]
    return np.linspace(scan["eta_min"], scan["eta_max"], scan["eta_points"])


def optimum_window(curve: pd.DataFrame, fraction: float = 0.95) -> dict[str, float]:
    keep = curve[curve.zt >= fraction * curve.zt.max()]
    best = curve.loc[curve.zt.idxmax()]
    result = {"zt_max": best.zt, "eta_opt": best.eta, "fraction": fraction}
    for col in (
        "n_sheet_cm2",
        "n_equiv_cm3",
        "seebeck_uV_K",
        "sigma_S_m",
        "lorenz_W_ohm_K2",
        "kappa_e_W_mK",
    ):
        result[f"{col}_opt"] = best[col]
        result[f"{col}_low"] = keep[col].min()
        result[f"{col}_high"] = keep[col].max()
    return result


def build_baseline(cfg, ep, lp, sp, eta):
    lattice = lattice_transport(lp, sp)
    curve = pd.DataFrame(electronic_transport(eta, ep, lattice["kappa_lattice_W_mK"], sp))
    curve.to_csv(OUT / "baseline_eta_scan.csv", index=False)
    window = optimum_window(curve)
    pd.DataFrame([window]).to_csv(OUT / "baseline_optimum_window.csv", index=False)

    rows = []
    for temperature in (300.0, 450.0, 600.0, 750.0, 900.0):
        ep_t = replace(ep, temperature_K=temperature)
        lp_t = replace(lp, temperature_K=temperature)
        opt = optimise_eta(eta, ep_t, lp_t, sp)
        opt.update(
            {
                "temperature_K": temperature,
                "gap_8kBT_eV": bipolar_gap_threshold_eV(temperature, 8),
                "gap_10kBT_eV": bipolar_gap_threshold_eV(temperature, 10),
            }
        )
        rows.append(opt)
    temperature = pd.DataFrame(rows)
    temperature.to_csv(OUT / "temperature_optima.csv", index=False)
    return curve, window, temperature, lattice


def build_mass_surfaces(ep, lp, sp, eta):
    md_values = np.linspace(0.2, 2.5, 42)
    mc_values = np.linspace(0.10, 1.20, 42)
    zt = np.empty((len(mc_values), len(md_values)))
    n_sheet = np.empty_like(zt)
    eta_opt = np.empty_like(zt)
    for i, mc in enumerate(mc_values):
        for j, md in enumerate(md_values):
            opt = optimise_eta(eta, replace(ep, dos_mass_me=md, conductivity_mass_me=mc), lp, sp)
            zt[i, j] = opt["zt"]
            n_sheet[i, j] = opt["n_sheet_cm2"]
            eta_opt[i, j] = opt["eta"]
    rows = []
    for i, mc in enumerate(mc_values):
        for j, md in enumerate(md_values):
            rows.append(
                {
                    "dos_mass_me": md,
                    "conductivity_mass_me": mc,
                    "zt_opt": zt[i, j],
                    "n_sheet_opt_cm2": n_sheet[i, j],
                    "eta_opt": eta_opt[i, j],
                }
            )
    pd.DataFrame(rows).to_csv(OUT / "mass_response_surface.csv", index=False)
    return md_values, mc_values, zt, n_sheet


def build_structure_surfaces(cfg, ep, lp, sp, eta):
    p_values = np.linspace(0.0, cfg["scan"]["porosity_max"], 45)
    w_values = np.linspace(0.0, cfg["scan"]["wrinkle_slope_max"], 45)
    zt = np.empty((len(w_values), len(p_values)))
    kl = np.empty_like(zt)
    sigma = np.empty_like(zt)
    rows = []
    for i, w in enumerate(w_values):
        for j, porosity in enumerate(p_values):
            local = replace(sp, porosity=porosity, wrinkle_slope=w)
            lattice = lattice_transport(lp, local)
            opt = optimise_eta(eta, ep, lp, local)
            zt[i, j] = opt["zt"]
            kl[i, j] = lattice["kappa_lattice_W_mK"]
            sigma[i, j] = opt["sigma_S_m"]
            rows.append(
                {
                    "porosity": porosity,
                    "wrinkle_slope": w,
                    "zt_opt": opt["zt"],
                    "eta_opt": opt["eta"],
                    "n_sheet_opt_cm2": opt["n_sheet_cm2"],
                    "sigma_opt_S_m": opt["sigma_S_m"],
                    **lattice,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "porosity_wrinkle_surface.csv", index=False)

    stiffness_values = np.linspace(0.25, 1.25, 45)
    overlap_values = np.linspace(0.0, 1.0, 45)
    zt_soft = np.empty((len(overlap_values), len(stiffness_values)))
    kl_soft = np.empty_like(zt_soft)
    soft_rows = []
    for i, overlap in enumerate(overlap_values):
        for j, stiffness in enumerate(stiffness_values):
            local = replace(sp, stiffness_ratio=stiffness, acoustic_optical_overlap=overlap)
            lattice = lattice_transport(lp, local)
            opt = optimise_eta(eta, ep, lp, local)
            zt_soft[i, j] = opt["zt"]
            kl_soft[i, j] = lattice["kappa_lattice_W_mK"]
            soft_rows.append(
                {
                    "stiffness_ratio": stiffness,
                    "acoustic_optical_overlap": overlap,
                    "zt_opt": opt["zt"],
                    "sigma_opt_S_m": opt["sigma_S_m"],
                    **lattice,
                }
            )
    pd.DataFrame(soft_rows).to_csv(OUT / "softness_overlap_surface.csv", index=False)
    return (
        p_values,
        w_values,
        zt,
        kl,
        sigma,
        stiffness_values,
        overlap_values,
        zt_soft,
        kl_soft,
        frame,
    )


def quadratic_fit(frame: pd.DataFrame):
    p = frame.porosity.to_numpy()
    w = frame.wrinkle_slope.to_numpy()
    y = frame.zt_opt.to_numpy()
    x = np.column_stack([np.ones_like(p), p, w, p**2, p * w, w**2])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ beta
    r2 = 1.0 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    names = ["intercept", "porosity", "wrinkle_slope", "porosity_sq", "interaction", "wrinkle_slope_sq"]
    coef = pd.DataFrame({"term": names, "coefficient": beta})
    coef["r_squared"] = r2
    coef.to_csv(OUT / "quadratic_surrogate_coefficients.csv", index=False)
    out = frame[["porosity", "wrinkle_slope", "zt_opt"]].copy()
    out["zt_quadratic"] = pred
    out["residual"] = y - pred
    out.to_csv(OUT / "quadratic_surrogate_predictions.csv", index=False)
    return beta, r2, out


def robust_sampling(cfg, ep, lp, sp, eta):
    rng = np.random.default_rng(cfg["scan"]["random_seed"])
    rows = []
    for _ in range(cfg["scan"]["robust_samples"]):
        local_ep = replace(
            ep,
            dos_mass_me=rng.uniform(0.35, 1.8),
            conductivity_mass_me=rng.uniform(0.12, 0.9),
            valley_degeneracy=int(rng.integers(1, 7)),
            elastic_modulus_2d_N_per_m=rng.uniform(25, 110),
            deformation_potential_eV=rng.uniform(3.0, 10.0),
            effective_thickness_nm=rng.uniform(0.65, 1.25),
            mobility_retention=rng.uniform(0.45, 0.95),
        )
        local_sp = replace(
            sp,
            porosity=rng.uniform(0.0, 0.30),
            pore_spacing_nm=rng.uniform(5.0, 80.0),
            wrinkle_slope=rng.uniform(0.0, 0.8),
            wrinkle_wavelength_nm=rng.uniform(5.0, 80.0),
            stiffness_ratio=rng.uniform(0.45, 1.15),
            mass_ratio=rng.uniform(0.7, 1.8),
            acoustic_optical_overlap=rng.uniform(0.05, 0.90),
        )
        local_lp = replace(
            lp,
            heat_capacity_J_per_m3K=rng.uniform(1.1e6, 2.2e6),
            sound_velocity_m_per_s=rng.uniform(1500, 4500),
            intrinsic_lifetime_ps_at_600K=rng.uniform(0.15, 1.2),
        )
        lattice = lattice_transport(local_lp, local_sp)
        opt = optimise_eta(eta, local_ep, local_lp, local_sp)
        rows.append({**asdict(local_ep), **asdict(local_sp), **lattice, **opt})
    samples = pd.DataFrame(rows)
    samples.to_csv(OUT / "robust_scenario_samples.csv", index=False)

    output_cols = [
        "eta",
        "n_sheet_cm2",
        "n_equiv_cm3",
        "seebeck_uV_K",
        "mobility_cm2_Vs",
        "sigma_S_m",
        "lorenz_W_ohm_K2",
        "kappa_e_W_mK",
        "kappa_lattice_W_mK",
        "group_velocity_m_s",
        "lifetime_ps",
        "mean_free_path_nm",
        "zt",
    ]
    quantiles = samples[output_cols].quantile([0.1, 0.5, 0.9]).T.reset_index()
    quantiles.columns = ["quantity", "p10", "median", "p90"]
    quantiles.to_csv(OUT / "robust_optimum_quantiles.csv", index=False)

    top = samples[samples.zt >= samples.zt.quantile(0.9)]
    design_cols = [
        "dos_mass_me",
        "conductivity_mass_me",
        "valley_degeneracy",
        "elastic_modulus_2d_N_per_m",
        "deformation_potential_eV",
        "porosity",
        "wrinkle_slope",
        "stiffness_ratio",
        "mass_ratio",
        "acoustic_optical_overlap",
        "group_velocity_m_s",
        "lifetime_ps",
        "mean_free_path_nm",
        "kappa_lattice_W_mK",
    ]
    top_ranges = top[design_cols].quantile([0.1, 0.5, 0.9]).T.reset_index()
    top_ranges.columns = ["design_variable", "p10", "median", "p90"]
    top_ranges.to_csv(OUT / "top_decile_design_ranges.csv", index=False)
    return samples, quantiles, top_ranges


def jarvis_anchor():
    path = REPO / "jarvis_2d_te_atlas" / "data" / "processed" / "ZT_e_all.csv"
    if not path.exists():
        return pd.DataFrame()
    data = pd.read_csv(path)
    data["implied_mobility_cm2_Vs"] = (
        10 ** data.log_sigma_dom_geo / (1e26 * 1.602176634e-19) * 1e4
    )
    rows = []
    for carrier, group in data.groupby("carrier"):
        mass_col = "m_elec_median" if carrier == "n" else "m_hole_median"
        valid_mass = group.loc[(group[mass_col] > 0) & (group[mass_col] < 100), mass_col]
        rows.append(
            {
                "carrier": carrier,
                "count": len(group),
                "mobility_p25_cm2_Vs": group.implied_mobility_cm2_Vs.quantile(0.25),
                "mobility_median_cm2_Vs": group.implied_mobility_cm2_Vs.median(),
                "mobility_p75_cm2_Vs": group.implied_mobility_cm2_Vs.quantile(0.75),
                "effective_mass_median_me_filtered": valid_mass.median(),
                "note": "mobility inferred from JARVIS sigma at fixed n=1e20 cm^-3; not an independent measurement",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "jarvis_anchor_summary.csv", index=False)
    return out


def plot_formula_tree():
    fig, ax = plt.subplots(figsize=(12, 6.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")
    box = dict(boxstyle="round,pad=0.35", ec="#4a4a4a", lw=0.9)

    nodes = [
        (6, 6.35, r"$ZT=S^2\sigma T/(\kappa_e+\kappa_L)$", "#fff2b2"),
        (2.0, 5.0, r"$S(\eta,r,d)$", "#cfe8ff"),
        (4.5, 5.0, r"$\sigma=(n_{2D}/t)e\mu$", "#cfe8ff"),
        (7.2, 5.0, r"$\kappa_e=L(\eta,r)\sigma T$", "#cfe8ff"),
        (10.2, 5.0, r"$\kappa_L=\sum C_\lambda v_\lambda^2\tau_\lambda/V$", "#ffd8c9"),
        (1.6, 3.3, r"$\eta=(E_F-E_{edge})/k_BT$", "#e7f3ff"),
        (4.2, 3.3, r"$n_{2D}\propto N_vm_d^*TF_0(\eta)$", "#e7f3ff"),
        (6.8, 3.3, r"$\mu\propto C_{2D}/(Tm_c^*m_d^*E_1^2)$", "#e7f3ff"),
        (9.2, 3.3, r"$v_g=\partial\omega/\partial q$", "#ffe8df"),
        (11.0, 3.3, r"$\tau^{-1}=\sum_i\tau_i^{-1}$", "#ffe8df"),
        (2.3, 1.55, "doping / gate\ncarrier density", "#f4f4f4"),
        (5.2, 1.55, "valley count, DOS mass,\nconductivity mass, E1, C2D", "#f4f4f4"),
        (8.5, 1.55, "bond stiffness / mass\nacoustic-optical overlap", "#f4f4f4"),
        (11.0, 1.55, "pores / ligament scale\nwrinkles / boundaries", "#f4f4f4"),
    ]
    for x, y, text, color in nodes:
        ax.text(x, y, text, ha="center", va="center", bbox={**box, "fc": color}, fontsize=9.5)
    arrows = [
        ((6, 6.0), (2, 5.35)), ((6, 6.0), (4.5, 5.35)), ((6, 6.0), (7.2, 5.35)), ((6, 6.0), (10.2, 5.35)),
        ((2, 4.65), (1.6, 3.65)), ((4.5, 4.65), (4.2, 3.65)), ((4.5, 4.65), (6.8, 3.65)),
        ((10.2, 4.65), (9.2, 3.65)), ((10.2, 4.65), (11.0, 3.65)),
        ((1.6, 2.95), (2.3, 1.95)), ((4.2, 2.95), (5.2, 1.95)), ((6.8, 2.95), (5.2, 1.95)),
        ((9.2, 2.95), (8.5, 1.95)), ((11, 2.95), (11, 1.95)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#666", lw=1.0))
    ax.text(6, 0.45, "Exact identities  →  reduced transport variables  →  measurable structural proxies", ha="center", fontsize=10, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "01_formula_decomposition.png", bbox_inches="tight")
    plt.close(fig)


def plot_electronic(curve, md, mc, zt_mass, n_mass, temperature, ep):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    ax = axs[0, 0]
    ax.plot(curve.n_equiv_cm3, curve.zt, color="#c43d3d", lw=2, label="ZT")
    ax.set_xscale("log")
    ax.set_xlabel(r"equivalent carrier concentration $n_{eq}$ (cm$^{-3}$)")
    ax.set_ylabel("ZT", color="#c43d3d")
    ax.tick_params(axis="y", labelcolor="#c43d3d")
    bx = ax.twinx()
    bx.plot(curve.n_equiv_cm3, curve.seebeck_uV_K, color="#2867b2", lw=1.7, label="|S|")
    bx.set_ylabel(r"$|S|$ ($\mu$V/K)", color="#2867b2")
    bx.tick_params(axis="y", labelcolor="#2867b2")
    best = curve.loc[curve.zt.idxmax()]
    ax.axvline(best.n_equiv_cm3, color="k", ls="--", lw=0.9)
    ax.set_title("(a) Doping trade-off with variable Lorenz number")

    im = axs[0, 1].imshow(
        zt_mass,
        origin="lower",
        aspect="auto",
        extent=[md.min(), md.max(), mc.min(), mc.max()],
        cmap="viridis",
    )
    axs[0, 1].scatter([ep.dos_mass_me], [ep.conductivity_mass_me], marker="*", s=80, c="white", ec="k")
    axs[0, 1].set_xlabel(r"DOS mass $m_d^*/m_e$")
    axs[0, 1].set_ylabel(r"conductivity mass $m_c^*/m_e$")
    axs[0, 1].set_title("(b) Optimized ZT: low conductivity mass is the lever")
    fig.colorbar(im, ax=axs[0, 1], label="max ZT over doping")

    md_line = np.linspace(0.2, 2.5, 70)
    for nv, color in [(1, "#777777"), (2, "#2867b2"), (4, "#d17a00")]:
        values = []
        for md0 in md_line:
            p = replace(ep, dos_mass_me=md0, valley_degeneracy=nv)
            opt = optimise_eta(curve.eta.to_numpy(), p, LatticeParams(temperature_K=ep.temperature_K), StructureParams())
            values.append(opt["n_sheet_cm2"])
        axs[1, 0].plot(md_line, values, lw=1.8, color=color, label=f"Nv={nv}")
    axs[1, 0].set_yscale("log")
    axs[1, 0].set_xlabel(r"DOS mass $m_d^*/m_e$")
    axs[1, 0].set_ylabel(r"optimal sheet density $n_{2D}$ (cm$^{-2}$)")
    axs[1, 0].legend(frameon=False)
    axs[1, 0].set_title("(c) DOS mass and valley count set the required doping")

    s, lorenz = spb_seebeck_lorenz(curve.eta.to_numpy())
    axs[1, 1].plot(np.asarray(s) * 1e6, np.asarray(lorenz) * 1e8, color="#5d4a9c", lw=2)
    axs[1, 1].axhline(2.44, color="grey", ls="--", label=r"metal value $L_0$")
    axs[1, 1].set_xlabel(r"$|S|$ ($\mu$V/K)")
    axs[1, 1].set_ylabel(r"Lorenz number ($10^{-8}$ W$\Omega$K$^{-2}$)")
    axs[1, 1].legend(frameon=False)
    axs[1, 1].set_title("(d) Fixed Lorenz number is inaccurate near optimum")
    fig.suptitle("2D single-parabolic-band electronic design maps", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "02_electronic_design_maps.png", bbox_inches="tight")
    plt.close(fig)


def plot_lattice_structure(
    p_values,
    w_values,
    zt,
    stiffness,
    overlap,
    zt_soft,
    temperature,
    ep,
):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    velocities = np.linspace(1000, 5000, 100)
    lifetimes = np.geomspace(0.05, 2.0, 100)
    vv, tt = np.meshgrid(velocities, lifetimes)
    kk = 0.5 * 1.6e6 * vv**2 * tt * 1e-12
    im0 = axs[0, 0].contourf(velocities, lifetimes, kk, levels=np.linspace(0, 12, 25), cmap="magma_r", extend="max")
    axs[0, 0].set_yscale("log")
    axs[0, 0].contour(velocities, lifetimes, kk, levels=[0.5, 1.0, 1.5, 3.0], colors="white", linewidths=0.8)
    axs[0, 0].set_xlabel(r"projected group velocity $v_g$ (m/s)")
    axs[0, 0].set_ylabel(r"effective lifetime $\tau$ (ps)")
    axs[0, 0].set_title(r"(a) $\kappa_L=C_Vv_g^2\tau/2$")
    fig.colorbar(im0, ax=axs[0, 0], label=r"$\kappa_L$ (W m$^{-1}$ K$^{-1}$)")

    im1 = axs[0, 1].imshow(
        zt,
        origin="lower",
        aspect="auto",
        extent=[p_values.min(), p_values.max(), w_values.min(), w_values.max()],
        cmap="viridis",
    )
    idx = np.unravel_index(np.argmax(zt), zt.shape)
    axs[0, 1].scatter(p_values[idx[1]], w_values[idx[0]], marker="*", s=90, c="white", ec="k")
    axs[0, 1].set_xlabel("porosity")
    axs[0, 1].set_ylabel(r"wrinkle slope $2\pi A/\lambda$")
    axs[0, 1].set_title("(b) Coupled electronic/phonon response")
    fig.colorbar(im1, ax=axs[0, 1], label="max ZT over doping")

    im2 = axs[1, 0].imshow(
        zt_soft,
        origin="lower",
        aspect="auto",
        extent=[stiffness.min(), stiffness.max(), overlap.min(), overlap.max()],
        cmap="viridis",
    )
    idx2 = np.unravel_index(np.argmax(zt_soft), zt_soft.shape)
    axs[1, 0].scatter(stiffness[idx2[1]], overlap[idx2[0]], marker="*", s=90, c="white", ec="k")
    axs[1, 0].set_xlabel(r"stiffness ratio $C_{2D}/C_{ref}$")
    axs[1, 0].set_ylabel("acoustic-optical overlap index")
    axs[1, 0].set_title("(c) Softness is a trade-off; branch overlap targets phonons")
    fig.colorbar(im2, ax=axs[1, 0], label="max ZT over doping")

    ax = axs[1, 1]
    ax.plot(temperature.temperature_K, temperature.n_sheet_cm2, "-o", color="#2867b2", label=r"$n_{2D}^*$")
    ax.set_yscale("log")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"optimal $n_{2D}$ (cm$^{-2}$)", color="#2867b2")
    ax.tick_params(axis="y", labelcolor="#2867b2")
    bx = ax.twinx()
    bx.plot(temperature.temperature_K, temperature.gap_10kBT_eV, "--s", color="#c43d3d", label=r"$10k_BT$")
    bx.set_ylabel("conservative minimum gap (eV)", color="#c43d3d")
    bx.tick_params(axis="y", labelcolor="#c43d3d")
    ax.set_title("(d) Temperature shifts doping and bipolar-gap constraints")
    fig.suptitle("Lattice and structural response maps", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "03_lattice_structure_maps.png", bbox_inches="tight")
    plt.close(fig)


def plot_quadratic(p_values, w_values, zt, quadratic, beta, r2):
    fit = quadratic.zt_quadratic.to_numpy().reshape(len(w_values), len(p_values))
    residual = quadratic.residual.to_numpy().reshape(len(w_values), len(p_values))
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.5))
    extent = [p_values.min(), p_values.max(), w_values.min(), w_values.max()]
    im0 = axs[0, 0].imshow(zt, origin="lower", aspect="auto", extent=extent, cmap="viridis")
    axs[0, 0].set_title("(a) Reduced-physics response")
    fig.colorbar(im0, ax=axs[0, 0], label="ZT")
    im1 = axs[0, 1].imshow(fit, origin="lower", aspect="auto", extent=extent, cmap="viridis", vmin=zt.min(), vmax=zt.max())
    axs[0, 1].set_title(f"(b) Second-order surface, R²={r2:.4f}")
    fig.colorbar(im1, ax=axs[0, 1], label="ZT")
    lim = np.max(np.abs(residual))
    im2 = axs[1, 0].imshow(
        residual,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
    )
    axs[1, 0].set_title("(c) Physics model − quadratic surrogate")
    fig.colorbar(im2, ax=axs[1, 0], label="ZT residual")
    for ax in axs.ravel()[:3]:
        ax.set_xlabel("porosity")
        ax.set_ylabel("wrinkle slope")

    for target_w, color in [(0.0, "#2867b2"), (0.6, "#d17a00"), (1.2, "#c43d3d")]:
        i = np.argmin(np.abs(w_values - target_w))
        axs[1, 1].plot(p_values, zt[i], color=color, lw=2, label=f"physics, w={w_values[i]:.1f}")
        axs[1, 1].plot(p_values, fit[i], color=color, lw=1.2, ls="--")
    axs[1, 1].set_xlabel("porosity")
    axs[1, 1].set_ylabel("ZT")
    axs[1, 1].legend(frameon=False, fontsize=8)
    axs[1, 1].set_title("(d) Solid=physics, dashed=quadratic")
    fig.suptitle("Second-order response image: useful locally, not a new physical law", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "04_quadratic_response_surface.png", bbox_inches="tight")
    plt.close(fig)


def plot_feature_matrix(top_ranges):
    labels = [
        "wrinkle slope",
        "porosity",
        "small ligament",
        "softness",
        "heavy atoms",
        "A-O overlap",
        "valley degeneracy",
        "low conductivity mass",
    ]
    columns = [r"$n/S$", r"$\mu/\sigma$", r"$v_g$", r"$\tau$", r"$\kappa_L$", "ZT"]
    # +1 increase, -1 decrease, +/-0.5 conditional or indirect.
    matrix = np.array(
        [
            [0.0, -0.5, -1.0, -1.0, -1.0, 0.5],
            [0.0, -1.0, 0.0, -1.0, -1.0, 0.5],
            [0.0, -0.5, 0.0, -1.0, -1.0, 0.5],
            [0.0, -1.0, -1.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0, -0.5, -1.0, 0.5],
            [0.0, 0.0, -0.5, -1.0, -1.0, 1.0],
            [1.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4), gridspec_kw={"width_ratios": [1.1, 1]})
    im = axs[0].imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    axs[0].set_xticks(np.arange(len(columns)), columns)
    axs[0].set_yticks(np.arange(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = "+" if value == 1 else "−" if value == -1 else "cond." if abs(value) == 0.5 else "—"
            axs[0].text(j, i, text, ha="center", va="center", fontsize=8)
    axs[0].set_title("(a) Directional mechanism matrix")
    fig.colorbar(im, ax=axs[0], fraction=0.045, label="direction (conditional at ±0.5)")

    names = [
        "conductivity_mass_me",
        "deformation_potential_eV",
        "porosity",
        "wrinkle_slope",
        "stiffness_ratio",
        "acoustic_optical_overlap",
        "group_velocity_m_s",
        "mean_free_path_nm",
        "kappa_lattice_W_mK",
    ]
    subset = top_ranges.set_index("design_variable").loc[names].copy()
    y = np.arange(len(subset))
    width = subset.p90 - subset.p10
    axs[1].barh(y, width, left=subset.p10, color="#6aa6d8", alpha=0.75)
    axs[1].plot(subset["median"], y, "ko", ms=4)
    axs[1].set_yticks(y, names)
    axs[1].invert_yaxis()
    axs[1].set_xscale("symlog", linthresh=0.2)
    axs[1].set_xlabel("physical value (mixed units; log-like display)")
    axs[1].set_title("(b) Top-decile scenario ranges (10–90%, dot=median)")
    fig.tight_layout()
    fig.savefig(FIG / "05_feature_relationships_and_ranges.png", bbox_inches="tight")
    plt.close(fig)


def fmt(value, digits=3):
    if abs(value) >= 1e4 or (0 < abs(value) < 1e-3):
        return f"{value:.2e}"
    return f"{value:.{digits}g}"


def write_summary(
    ep,
    lattice,
    window,
    temperature,
    structure_frame,
    zt_soft,
    stiffness,
    overlap,
    beta,
    r2,
    quantiles,
    top_ranges,
    anchor,
):
    best_struct = structure_frame.loc[structure_frame.zt_opt.idxmax()]
    idx_soft = np.unravel_index(np.argmax(zt_soft), zt_soft.shape)
    q = quantiles.set_index("quantity")
    t = top_ranges.set_index("design_variable")
    lines = [
        "# 数值推演结果与设计窗口",
        "",
        "## 先给结论",
        "",
        f"基准情景（{ep.temperature_K:.0f} K）在完整掺杂扫描中的最高 `ZT={window['zt_max']:.3f}`；95% 峰值平台对应 "
        f"`n2D={fmt(window['n_sheet_cm2_low'])}–{fmt(window['n_sheet_cm2_high'])} cm^-2`，"
        f"体等效浓度 `neq={fmt(window['n_equiv_cm3_low'])}–{fmt(window['n_equiv_cm3_high'])} cm^-3`，"
        f"`|S|={window['seebeck_uV_K_low']:.0f}–{window['seebeck_uV_K_high']:.0f} μV/K`。",
        "",
        f"基准晶格侧给出 `κL={lattice['kappa_lattice_W_mK']:.3f} W m^-1 K^-1`、"
        f"投影群速度 `{lattice['group_velocity_m_s']:.0f} m/s`、"
        f"寿命 `{lattice['lifetime_ps']:.3f} ps` 和平均自由程 `{lattice['mean_free_path_nm']:.3f} nm`。"
        "这些数值是压缩模型的情景值，不是对某个材料的第一性原理预测。",
        "",
        "最重要的结论不是一个孤立的‘最佳有效质量’：二维形变势模型中，`md*` 提高会增加态密度和所需载流子数，但同时降低迁移率；在给定散射假设下两者大体抵消。真正有利的是较低的导电质量 `mc*`、较高的谷简并度 `Nv`，以及较小的形变势 `E1`。",
        "",
        "## 基准最优点",
        "",
        "| 量 | 最优点 | 95% 峰值平台 |",
        "|---|---:|---:|",
    ]
    labels = [
        ("约化费米能级 η", "eta", ""),
        ("n2D (cm^-2)", "n_sheet_cm2", ""),
        ("neq (cm^-3)", "n_equiv_cm3", ""),
        ("|S| (μV/K)", "seebeck_uV_K", ""),
        ("σ (S/m)", "sigma_S_m", ""),
        ("Lorenz (W Ω K^-2)", "lorenz_W_ohm_K2", ""),
        ("κe (W/mK)", "kappa_e_W_mK", ""),
    ]
    for label, key, _ in labels:
        if key == "eta":
            lines.append(f"| {label} | {fmt(window['eta_opt'])} | — |")
        else:
            lines.append(
                f"| {label} | {fmt(window[key + '_opt'])} | {fmt(window[key + '_low'])} – {fmt(window[key + '_high'])} |"
            )
    lines += [
        "",
        "## 跨情景稳健范围",
        "",
        "下表是 400 个透明参数情景分别优化掺杂后得到的 10–90% 区间，表示模型敏感性，不是统计置信区间：",
        "",
        "| 量 | P10 | 中位 | P90 |",
        "|---|---:|---:|---:|",
    ]
    for name in [
        "n_sheet_cm2",
        "n_equiv_cm3",
        "seebeck_uV_K",
        "mobility_cm2_Vs",
        "sigma_S_m",
        "lorenz_W_ohm_K2",
        "kappa_e_W_mK",
        "kappa_lattice_W_mK",
        "group_velocity_m_s",
        "lifetime_ps",
        "mean_free_path_nm",
        "zt",
    ]:
        row = q.loc[name]
        lines.append(f"| {name} | {fmt(row.p10)} | {fmt(row['median'])} | {fmt(row.p90)} |")
    lines += [
        "",
        "## 高表现情景的条件范围",
        "",
        "按每个情景都已优化掺杂后的 ZT 取前 10%，其输入/中间量 10–90% 范围如下。由于采样不是材料数据库，边界值只用于设定筛选盒：",
        "",
        "| 设计量 | P10 | 中位 | P90 |",
        "|---|---:|---:|---:|",
    ]
    for name, row in t.iterrows():
        lines.append(f"| {name} | {fmt(row.p10)} | {fmt(row['median'])} | {fmt(row.p90)} |")
    lines += [
        "",
        "## 结构响应与二次图像",
        "",
        f"在基准情景的孔隙率—褶皱斜率网格内，网格最高点为 `φ={best_struct.porosity:.3f}`、"
        f"`2πA/λ={best_struct.wrinkle_slope:.3f}`、`ZT={best_struct.zt_opt:.3f}`。"
        "若最高点落在扫描边界，它只能说明当前假设下仍在单调变化，不能宣称存在内部最优。",
        "",
        f"刚度—声光重叠网格的最高点为 `stiffness_ratio={stiffness[idx_soft[1]]:.3f}`、"
        f"`overlap={overlap[idx_soft[0]]:.3f}`、`ZT={zt_soft[idx_soft]:.3f}`。"
        "软化会同时降低群速度和形变势迁移率，因此不是越软越好；提高声—光支耦合在本模型中主要缩短声子寿命，对电子侧没有直接奖励。",
        "",
        f"孔隙率—褶皱面用二次多项式拟合得到 `R²={r2:.5f}`。系数保存在 `quadratic_surrogate_coefficients.csv`；"
        "二次面只是扫描域内的局部代理，不是普适物理定律。",
        "",
        "## 温度与带隙",
        "",
        "| T (K) | n2D* (cm^-2) | η* | 保守 Eg,min=10kBT (eV) |",
        "|---:|---:|---:|---:|",
    ]
    for _, row in temperature.iterrows():
        lines.append(f"| {row.temperature_K:.0f} | {fmt(row.n_sheet_cm2)} | {row.eta:.2f} | {row.gap_10kBT_eV:.3f} |")
    lines += [
        "",
        "SPB 只能给出抑制双极输运所需的最小带隙尺度，不能推出通用的 1–2 eV 上下界。过大的带隙是否不利取决于它与质量、谷简并度、缺陷能级和迁移率的材料相关共变。",
        "",
        "## 特征取舍",
        "",
        "- 保留并量化：褶皱斜率/波长、孔隙率/骨架尺度、刚度比、质量比、声—光支重叠或避免交叉强度。",
        "- 替换：‘骨架’拆成孔隙率、连通性、迂曲度和韧带尺度；‘振动频率’拆成群速度、模态热容和寿命。",
        "- 删除：‘结构是否均质’二元标签、声学/光学支的简单数量。它们都没有唯一的单调物理映射。",
        "",
        "## 与原项目数据的锚定",
        "",
    ]
    if not anchor.empty:
        lines += [
            "JARVIS 输运表在固定 `n=1e20 cm^-3, T=600 K` 条件下可反推一个等效迁移率，仅用于校准数量级：",
            "",
            "| carrier | μ P25 | μ median | μ P75 | 过滤后质量中位 (me) |",
            "|---|---:|---:|---:|---:|",
        ]
        for _, row in anchor.iterrows():
            lines.append(
                f"| {row.carrier} | {row.mobility_p25_cm2_Vs:.2f} | {row.mobility_median_cm2_Vs:.2f} | "
                f"{row.mobility_p75_cm2_Vs:.2f} | {row.effective_mass_median_me_filtered:.2f} |"
            )
        lines.append("")
        lines.append("该迁移率是由 σ/(ne) 反推，不是独立 Hall 测量；原始表的极端值也说明必须使用分位数而非均值。")
    lines += [
        "",
        "## 图表",
        "",
        "1. `figures/01_formula_decomposition.png`：公式到结构代理的完整链条；",
        "2. `figures/02_electronic_design_maps.png`：掺杂、两类质量、谷简并和 Lorenz 数；",
        "3. `figures/03_lattice_structure_maps.png`：群速度—寿命、孔隙—褶皱、软化—声光重叠；",
        "4. `figures/04_quadratic_response_surface.png`：真实低阶模型与二次代理及残差；",
        "5. `figures/05_feature_relationships_and_ranges.png`：结构特征方向矩阵与高表现范围。",
        "",
        "## 下一层计算条件",
        "",
        "要把这里的范围真正落实到一个候选结构，应按以下顺序补数据：",
        "",
        "1. 电子结构：`mx,my`、谷数、带隙和应变下带边位移 `E1`；",
        "2. 弹性：方向相关 `C2D`；",
        "3. 声子二阶力常数：色散、虚频、群速度、声—光支重叠；",
        "4. 声子三阶力常数：模态寿命和 κL；",
        "5. 对褶皱/孔洞超胞同时复算电子迁移率保持率，不能只计算 κL 的收益。",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    cfg, ep, lp, sp = load_params()
    eta = eta_grid(cfg)
    curve, window, temperature, lattice = build_baseline(cfg, ep, lp, sp, eta)
    md, mc, zt_mass, n_mass = build_mass_surfaces(ep, lp, sp, eta)
    (
        p_values,
        w_values,
        zt_struct,
        kl_struct,
        sigma_struct,
        stiffness,
        overlap,
        zt_soft,
        kl_soft,
        structure_frame,
    ) = build_structure_surfaces(cfg, ep, lp, sp, eta)
    beta, r2, quadratic = quadratic_fit(structure_frame)
    samples, quantiles, top_ranges = robust_sampling(cfg, ep, lp, sp, eta)
    anchor = jarvis_anchor()

    plot_formula_tree()
    plot_electronic(curve, md, mc, zt_mass, n_mass, temperature, ep)
    plot_lattice_structure(p_values, w_values, zt_struct, stiffness, overlap, zt_soft, temperature, ep)
    plot_quadratic(p_values, w_values, zt_struct, quadratic, beta, r2)
    plot_feature_matrix(top_ranges)
    write_summary(
        ep,
        lattice,
        window,
        temperature,
        structure_frame,
        zt_soft,
        stiffness,
        overlap,
        beta,
        r2,
        quantiles,
        top_ranges,
        anchor,
    )
    print(f"written: {OUT / 'summary.md'}")
    print(f"baseline max ZT: {window['zt_max']:.4f}")
    print(f"quadratic R2: {r2:.6f}")


if __name__ == "__main__":
    main()

