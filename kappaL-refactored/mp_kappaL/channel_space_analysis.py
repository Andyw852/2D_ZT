"""Validate property-space relations to electronic and lattice conductivity.

This module deliberately separates three questions:

1. Direct monotone association (Spearman, formula bootstrap confidence interval).
2. Global geometry (correlation between pairwise distances in a property space and
   distances in log conductivity).
3. Local geometry (k-nearest-neighbour overlap, with the common-cohort null).

The JARVIS electronic fields are constant-relaxation-time quantities at 600 K.
The MP Snyder field is a 300 K analytic model proxy.  Their formula-level join is
therefore diagnostic rather than a material-resolved thermoelectric measurement.
Experimental lattice conductivity is analysed separately and never copied onto
MP polymorphs.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

from mp_kappaL.metrics import (
    benjamini_hochberg,
    crossview_overlap,
    distance_spearman,
    knn_neighbor_matrix,
    pairwise_distance_matrix,
)


ROOT = Path(__file__).resolve().parent
PROC = ROOT / "processed"
FIG = ROOT / "figures"
SEED = 20260828
K = 10
MAX_GEOMETRY_N = 2500


def _seed(label: str) -> int:
    return int.from_bytes(
        hashlib.blake2b(label.encode("utf-8"), digest_size=4).digest(), "little"
    )


def _stable_sample(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    """Order-independent deterministic formula sample."""
    if len(frame) <= n:
        return frame.sort_values("canon").reset_index(drop=True)
    score = frame["canon"].map(
        lambda x: int.from_bytes(
            hashlib.blake2b(str(x).encode("utf-8"), digest_size=8).digest(), "little"
        )
    )
    return frame.loc[score.nsmallest(n).index].sort_values("canon").reset_index(drop=True)


def _finite(frame: pd.DataFrame, columns: list[str], positive: list[str] | None = None):
    mask = np.ones(len(frame), dtype=bool)
    for col in columns:
        mask &= np.isfinite(pd.to_numeric(frame[col], errors="coerce"))
    for col in positive or []:
        mask &= pd.to_numeric(frame[col], errors="coerce") > 0
    return frame.loc[mask].copy()


def _bootstrap_spearman(x, y, label: str, n_boot: int = 500):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(SEED + _seed(label))
    boots = np.empty(n_boot)
    for b in range(n_boot):
        ix = rng.integers(0, len(x), len(x))
        boots[b] = stats.spearmanr(x[ix], y[ix]).statistic
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(rho), float(p), float(lo), float(hi)


def _property_values(frame: pd.DataFrame, carrier: str):
    """Property label -> (values, provenance class)."""
    return {
        "band gap": (frame["Eg_opt"], "same-source electronic"),
        "|Seebeck|": (frame[f"S_{carrier}"].abs(), "same transport calculation"),
        "log conductivity": (
            np.log10(frame[f"sigma_{carrier}"].clip(lower=1e-30)),
            "same transport calculation; WF-related",
        ),
        "log power factor": (
            np.log10(frame[f"PF_{carrier}"].clip(lower=1e-30)),
            "derived from Seebeck and conductivity",
        ),
        "log bulk modulus": (
            np.log10(frame["bulk_vrh"].clip(lower=1e-30)),
            "MP lattice descriptor",
        ),
        "log shear modulus": (
            np.log10(frame["shear_vrh"].clip(lower=1e-30)),
            "MP lattice descriptor",
        ),
        "log Debye temperature": (
            np.log10(frame["debye"].clip(lower=1e-30)),
            "MP lattice descriptor; Snyder input",
        ),
        "log density": (
            np.log10(frame["density"].clip(lower=1e-30)),
            "MP lattice descriptor; Snyder-related",
        ),
        "formation energy": (frame["fe"], "cross-source chemistry proxy"),
    }


def direct_correlations(base: pd.DataFrame, experimental: pd.DataFrame):
    """Property-to-target correlations, preserving the target's provenance."""
    rows = []
    exp_join = base.merge(
        experimental[["canon", "kappa_L", "n_curves"]], on="canon", how="inner"
    )
    for carrier in ("n", "p"):
        ke = f"kappa_e_{carrier}"
        targets = [
            ("electronic kappa (CRTA)", base, ke, "JARVIS 600 K CRTA"),
            ("Snyder lattice kappa", base, "kappa_L_snyder_median",
             "MP 300 K model proxy"),
            ("experimental lattice kappa", exp_join, "kappa_L",
             "Starrydata formula-level experiment"),
        ]
        positive_property_source = {
            "log conductivity": f"sigma_{carrier}",
            "log power factor": f"PF_{carrier}",
            "log bulk modulus": "bulk_vrh",
            "log shear modulus": "shear_vrh",
            "log Debye temperature": "debye",
            "log density": "density",
        }
        for target, frame, target_column, target_source in targets:
            target_raw = pd.to_numeric(frame[target_column], errors="coerce")
            for prop, (x_all, relation) in _property_values(frame, carrier).items():
                valid = np.isfinite(pd.to_numeric(x_all, errors="coerce"))
                valid &= np.isfinite(target_raw) & (target_raw > 0)
                source_column = positive_property_source.get(prop)
                if source_column is not None:
                    valid &= pd.to_numeric(frame[source_column], errors="coerce") > 0
                x = np.asarray(x_all[valid], dtype=float)
                y = np.log10(target_raw[valid].to_numpy(float))
                rho, p, lo, hi = _bootstrap_spearman(
                    x, y, f"{carrier}-{target}-{prop}"
                )
                rows.append({
                    "carrier": carrier,
                    "target": target,
                    "property": prop,
                    "n": int(valid.sum()),
                    "spearman": rho,
                    "p_value": p,
                    "ci_lo": lo,
                    "ci_hi": hi,
                    "property_relation": relation,
                    "target_source": target_source,
                })
    out = pd.DataFrame(rows)
    out["q_bh_within_target"] = np.nan
    for _, ix in out.groupby(["carrier", "target"]).groups.items():
        out.loc[ix, "q_bh_within_target"] = benjamini_hochberg(
            out.loc[ix, "p_value"].to_numpy()
        )
    out.to_csv(PROC / "channel_property_correlations.csv", index=False)

    cross = []
    for carrier in ("n", "p"):
        ke = f"kappa_e_{carrier}"
        for target, cohort, kcol, source in [
            ("Snyder lattice kappa", base, "kappa_L_snyder_median", "model proxy"),
            ("experimental lattice kappa", exp_join, "kappa_L", "experiment"),
        ]:
            d = _finite(cohort, [ke, kcol], [ke, kcol])
            rho, p, lo, hi = _bootstrap_spearman(
                np.log10(d[ke]), np.log10(d[kcol]), f"cross-{carrier}-{target}"
            )
            cross.append({
                "carrier": carrier,
                "lattice_target": target,
                "source": source,
                "n": len(d),
                "spearman": rho,
                "p_value": p,
                "ci_lo": lo,
                "ci_hi": hi,
                "join_level": "canonical formula",
                "same_temperature": False,
                "material_id_resolved": False,
            })
    cross = pd.DataFrame(cross)
    cross["q_bh"] = benjamini_hochberg(cross["p_value"])
    cross.to_csv(PROC / "ke_kl_crosschannel.csv", index=False)
    return out, cross, exp_join


def experimental_subset_sensitivity(exp_join: pd.DataFrame):
    """Show how the small formula-level experimental result changes by match quality."""
    subsets = {
        "all formula matches": np.ones(len(exp_join), dtype=bool),
        "one MP polymorph": exp_join["n_mp_polymorphs"].eq(1).to_numpy(),
        "one experimental curve": exp_join["n_curves"].eq(1).to_numpy(),
        "one polymorph and one curve": (
            exp_join["n_mp_polymorphs"].eq(1) & exp_join["n_curves"].eq(1)
        ).to_numpy(),
    }
    properties = {
        "n-type electronic kappa": "kappa_e_n",
        "p-type electronic kappa": "kappa_e_p",
        "Debye temperature": "debye",
        "shear modulus": "shear_vrh",
        "bulk modulus": "bulk_vrh",
    }
    rows = []
    for subset_name, selection in subsets.items():
        subset = exp_join.loc[selection]
        for label, column in properties.items():
            cohort = _finite(subset, [column, "kappa_L"], [column, "kappa_L"])
            rho, p, lo, hi = _bootstrap_spearman(
                np.log10(cohort[column]), np.log10(cohort["kappa_L"]),
                f"subset-{subset_name}-{label}", n_boot=1000,
            )
            jackknife = []
            if len(cohort) >= 4:
                x = np.log10(cohort[column].to_numpy(float))
                y = np.log10(cohort["kappa_L"].to_numpy(float))
                for i in range(len(cohort)):
                    keep = np.arange(len(cohort)) != i
                    jackknife.append(stats.spearmanr(x[keep], y[keep]).statistic)
            rows.append({
                "subset": subset_name, "property": label, "n": len(cohort),
                "spearman": rho, "p_value": p, "ci_lo": lo, "ci_hi": hi,
                "leave_one_out_min": np.nanmin(jackknife),
                "leave_one_out_max": np.nanmax(jackknife),
            })
    out = pd.DataFrame(rows)
    out.to_csv(PROC / "experimental_subset_sensitivity.csv", index=False)
    return out


def _space_arrays(frame: pd.DataFrame, carrier: str):
    s = frame[f"S_{carrier}"].abs().to_numpy(float)
    eg = frame["Eg_opt"].to_numpy(float)
    logsigma = np.log10(frame[f"sigma_{carrier}"].to_numpy(float))
    lattice = np.column_stack([
        np.log10(frame["bulk_vrh"]),
        np.log10(frame["shear_vrh"]),
        np.log10(frame["debye"]),
        np.log10(frame["density"]),
    ])
    electronic_no_sigma = np.column_stack([eg, s])
    electronic = np.column_stack([eg, s, logsigma])
    return {
        "electronic without sigma": electronic_no_sigma,
        "electronic with sigma": electronic,
        "lattice descriptors": lattice,
        "joint descriptors": np.column_stack([electronic, lattice]),
    }


def _winsor_robust(x):
    x = np.asarray(x, dtype=float)
    lo = np.nanpercentile(x, 1, axis=0)
    hi = np.nanpercentile(x, 99, axis=0)
    return RobustScaler().fit_transform(np.clip(x, lo, hi))


def _heatmap(ax, frame, cmap, vmin, vmax, cbar_label, fmt=".2f"):
    """Small dependency-free annotated heatmap."""
    values = frame.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)
    image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(frame.shape[1]), labels=frame.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(frame.shape[0]), labels=frame.index)
    ax.set_xticks(np.arange(-.5, frame.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, frame.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    threshold = (vmin + vmax) / 2
    for i in range(frame.shape[0]):
        for j in range(frame.shape[1]):
            value = values[i, j]
            if np.isfinite(value):
                colour = "white" if value > threshold + .2 * (vmax - vmin) else "black"
                ax.text(j, i, format(value, fmt), ha="center", va="center",
                        fontsize=8.5, color=colour)
    cbar = ax.figure.colorbar(image, ax=ax, pad=.015)
    cbar.set_label(cbar_label)
    return image


def _geometry_one(cohort, carrier, target_name, target_values, provenance):
    ids = cohort["canon"].astype(str).to_numpy()
    target_distance = pairwise_distance_matrix(
        np.asarray(target_values, dtype=float).reshape(-1, 1), scale=False
    )
    target_nn = knn_neighbor_matrix(target_distance, K, ids=ids)
    drows, orows = [], []
    for space_name, values in _space_arrays(cohort, carrier).items():
        space_distance = pairwise_distance_matrix(_winsor_robust(values), scale=False)
        space_nn = knn_neighbor_matrix(space_distance, K, ids=ids)
        rng_d = np.random.RandomState(SEED + _seed(f"d-{carrier}-{target_name}-{space_name}"))
        dc = distance_spearman(
            space_distance, target_distance, rng_d,
            n_samp=50_000, n_mantel=30, n_boot=60,
        )
        drows.append({
            "carrier": carrier, "target": target_name, "space": space_name,
            "target_provenance": provenance, **dc,
        })
        rng_o = np.random.RandomState(SEED + _seed(f"o-{carrier}-{target_name}-{space_name}"))
        ov = crossview_overlap(space_nn, target_nn, K, rng_o, n_perm=300, n_boot=400)
        orows.append({
            "carrier": carrier, "target": target_name, "space": space_name,
            "target_provenance": provenance, **ov,
        })
    return drows, orows


def geometry_analysis(base: pd.DataFrame, exp_join: pd.DataFrame):
    drows, orows = [], []
    required_lattice = ["bulk_vrh", "shear_vrh", "debye", "density"]
    for carrier in ("n", "p"):
        ke = f"kappa_e_{carrier}"
        required = [ke, "kappa_L_snyder_median", "Eg_opt", f"S_{carrier}",
                    f"sigma_{carrier}"] + required_lattice
        cohort = _finite(
            base, required,
            [ke, "kappa_L_snyder_median", f"sigma_{carrier}"] + required_lattice,
        )
        cohort = _stable_sample(cohort, MAX_GEOMETRY_N)
        for target, y, source in [
            ("electronic kappa (CRTA)", np.log10(cohort[ke]), "JARVIS 600 K CRTA"),
            ("Snyder lattice kappa", np.log10(cohort["kappa_L_snyder_median"]),
             "MP 300 K model proxy"),
        ]:
            dr, ore = _geometry_one(cohort, carrier, target, y, source)
            drows.extend(dr)
            orows.extend(ore)

        exp_cohort = _finite(
            exp_join,
            [ke, "kappa_L", "Eg_opt", f"S_{carrier}", f"sigma_{carrier}"]
            + required_lattice,
            [ke, "kappa_L", f"sigma_{carrier}"] + required_lattice,
        )
        if len(exp_cohort) > K + 2:
            dr, ore = _geometry_one(
                exp_cohort.sort_values("canon").reset_index(drop=True), carrier,
                "experimental lattice kappa", np.log10(exp_cohort["kappa_L"]),
                "Starrydata formula-level experiment",
            )
            drows.extend(dr)
            orows.extend(ore)
    dframe = pd.DataFrame(drows)
    oframe = pd.DataFrame(orows)
    dframe.to_csv(PROC / "channel_space_distance.csv", index=False)
    oframe.to_csv(PROC / "channel_space_overlap.csv", index=False)
    return dframe, oframe


def plot_property_heatmap(correlations):
    labels = {
        "electronic kappa (CRTA)": "κe (CRTA)",
        "Snyder lattice kappa": "κL Snyder",
        "experimental lattice kappa": "κL experiment",
    }
    tmp = correlations.copy()
    tmp["column"] = tmp["carrier"] + " | " + tmp["target"].map(labels)
    order = [
        "band gap", "|Seebeck|", "log conductivity", "log power factor",
        "log bulk modulus", "log shear modulus", "log Debye temperature",
        "log density", "formation energy",
    ]
    cols = [f"{c} | {labels[t]}" for c in ("n", "p") for t in labels]
    matrix = tmp.pivot(index="property", columns="column", values="spearman").reindex(
        index=order, columns=cols
    )
    nmat = tmp.pivot(index="property", columns="column", values="n").reindex(
        index=order, columns=cols
    )
    fig, ax = plt.subplots(figsize=(13.5, 7.3))
    image = ax.imshow(matrix.to_numpy(float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(matrix.shape[1]), labels=matrix.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]), labels=matrix.index)
    ax.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i, row in enumerate(matrix.index):
        for j, col in enumerate(matrix.columns):
            value = matrix.loc[row, col]
            if np.isfinite(value):
                ax.text(j, i, f"{value:+.2f}\nN={int(nmat.loc[row, col])}",
                        ha="center", va="center", fontsize=7.5,
                        color="white" if abs(value) > .55 else "black")
    cbar = fig.colorbar(image, ax=ax, pad=.015)
    cbar.set_label("Spearman ρ")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Direct property–conductivity associations (formula level)")
    ax.text(
        0, -0.19,
        "κe: JARVIS 600 K CRTA. Snyder κL: MP 300 K analytic proxy. "
        "Experimental κL: Starrydata formula join; descriptor pairs N=96, electronic pairs N=70; "
        "polymorph-unresolved.",
        transform=ax.transAxes, fontsize=9, color="#444444",
    )
    fig.tight_layout()
    fig.savefig(FIG / "channel_property_correlations.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_geometry(dframe, oframe):
    d = dframe.copy()
    o = oframe.copy()
    d["column"] = d["carrier"] + " | " + d["target"].replace({
        "electronic kappa (CRTA)": "κe", "Snyder lattice kappa": "κL model",
        "experimental lattice kappa": "κL exp",
    })
    o["column"] = o["carrier"] + " | " + o["target"].replace({
        "electronic kappa (CRTA)": "κe", "Snyder lattice kappa": "κL model",
        "experimental lattice kappa": "κL exp",
    })
    spaces = ["electronic without sigma", "electronic with sigma",
              "lattice descriptors", "joint descriptors"]
    cols = [f"{c} | {t}" for c in ("n", "p") for t in ("κe", "κL model", "κL exp")]
    dm = d.pivot(index="space", columns="column", values="spearman").reindex(spaces, columns=cols)
    om = o.pivot(index="space", columns="column", values="enrichment").reindex(spaces, columns=cols)
    fig, axes = plt.subplots(2, 1, figsize=(13.2, 7.5), sharex=True)
    _heatmap(axes[0], dm, "viridis", 0, .85, "distance Spearman ρ")
    _heatmap(axes[1], om, "inferno_r", 0, max(3, np.nanpercentile(om, 95)),
             "kNN overlap / cohort null")
    axes[0].set_title("Global distance geometry")
    axes[1].set_title("Local 10-neighbour enrichment")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
    axes[1].text(
        0, -0.34,
        "Target-free property spaces; common cohort is recomputed before kNN. "
        "Full model cohorts N=2500; experimental cohorts N=70.",
        transform=axes[1].transAxes, fontsize=9, color="#444444",
    )
    fig.tight_layout()
    fig.savefig(FIG / "channel_space_geometry.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_pca(base):
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.4), constrained_layout=True)
    for row, carrier in enumerate(("n", "p")):
        ke = f"kappa_e_{carrier}"
        required = [ke, "kappa_L_snyder_median", "Eg_opt", f"S_{carrier}",
                    f"sigma_{carrier}", "bulk_vrh", "shear_vrh", "debye", "density"]
        cohort = _finite(
            base, required,
            [ke, "kappa_L_snyder_median", f"sigma_{carrier}", "bulk_vrh",
             "shear_vrh", "debye", "density"],
        )
        cohort = _stable_sample(cohort, 3000)
        joint = _space_arrays(cohort, carrier)["joint descriptors"]
        xy = PCA(n_components=2, random_state=SEED).fit_transform(_winsor_robust(joint))
        for col, (value, title, cmap) in enumerate([
            (np.log10(cohort[ke]), f"{carrier}-type | log κe (CRTA)", "viridis"),
            (np.log10(cohort["kappa_L_snyder_median"]),
             f"{carrier}-type | log κL (Snyder model)", "magma"),
        ]):
            ax = axes[row, col]
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=value, s=8, alpha=.72,
                            cmap=cmap, linewidths=0, rasterized=True)
            fig.colorbar(sc, ax=ax, pad=.01)
            ax.set_title(title)
            ax.set_xlabel("joint-property PC1")
            ax.set_ylabel("joint-property PC2")
    fig.suptitle("The same target-free joint property space, coloured by conductivity", fontsize=14)
    fig.savefig(FIG / "channel_space_pca.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_crosschannel(base, exp_join, cross):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.1), constrained_layout=True)
    for row, carrier in enumerate(("n", "p")):
        ke = f"kappa_e_{carrier}"
        full = _finite(base, [ke, "kappa_L_snyder_median"], [ke, "kappa_L_snyder_median"])
        exp = _finite(exp_join, [ke, "kappa_L"], [ke, "kappa_L"])
        r_model = cross[(cross.carrier == carrier) & (cross.source == "model proxy")].iloc[0]
        r_exp = cross[(cross.carrier == carrier) & (cross.source == "experiment")].iloc[0]
        ax = axes[row, 0]
        hb = ax.hexbin(np.log10(full[ke]), np.log10(full["kappa_L_snyder_median"]),
                       gridsize=45, mincnt=1, bins="log", cmap="Blues")
        fig.colorbar(hb, ax=ax, pad=.01, label="log count")
        ax.set_title(f"{carrier}-type: model proxy, ρ={r_model.spearman:+.2f} (N={r_model.n})")
        ax.set_xlabel("log κe (600 K CRTA scale)")
        ax.set_ylabel("log κL (300 K Snyder model)")

        ax = axes[row, 1]
        ax.scatter(np.log10(exp[ke]), np.log10(exp["kappa_L"]), s=24, alpha=.72,
                   color="#d95f02", edgecolor="white", linewidth=.3)
        ax.set_title(
            f"{carrier}-type: experiment, ρ={r_exp.spearman:+.2f} "
            f"[{r_exp.ci_lo:+.2f}, {r_exp.ci_hi:+.2f}] (N={r_exp.n})"
        )
        ax.set_xlabel("log κe (600 K CRTA scale)")
        ax.set_ylabel("log experimental κL (mixed T)")
    fig.suptitle("κe–κL cross-channel relation: proxy signal does not reproduce in experiment",
                 fontsize=14)
    fig.savefig(FIG / "ke_kl_crosschannel.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    base = pd.read_parquet(PROC / "dual_channel.parquet")
    experimental = pd.read_csv(PROC / "experimental_formula_targets.csv")
    correlations, cross, exp_join = direct_correlations(base, experimental)
    sensitivity = experimental_subset_sensitivity(exp_join)
    print("direct property correlations written")
    print(cross.to_string(index=False))
    print("\nexperimental subset sensitivity")
    print(sensitivity.to_string(index=False))
    distances, overlaps = geometry_analysis(base, exp_join)
    print("\nspace distance correlations")
    print(distances[["carrier", "target", "space", "n", "spearman", "ci_lo", "ci_hi"]].to_string(index=False))
    print("\nspace kNN overlap")
    print(overlaps[["carrier", "target", "space", "n", "overlap", "null_mean", "enrichment", "z"]].to_string(index=False))
    plot_property_heatmap(correlations)
    plot_geometry(distances, overlaps)
    plot_pca(base)
    plot_crosschannel(base, exp_join, cross)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
