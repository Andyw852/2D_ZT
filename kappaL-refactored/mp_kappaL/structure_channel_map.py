"""Map electronic- and lattice-favoured formulas in one target-free structure space.

The embedding uses composition-blind geometric SOAP only.  Conductivity targets are
used *after* PCA solely to colour points.  To avoid assigning one JARVIS formula to
several MP structures, the analysis is restricted to formulas having exactly one MP
structure in the cleaned elasticity cohort.

Two definitions are shown:

* literal: bottom-decile kappa_e versus bottom-decile Snyder-model kappa_L;
* thermoelectric: top-decile PF*T/kappa_e versus bottom-decile Snyder-model kappa_L.

The large-N lattice label remains a 300 K analytic model proxy.  Formula-matched
experimental low-kappa_L points are drawn as hollow stars for external orientation.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core import Composition
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parent
PROC = ROOT / "processed"
FIG = ROOT / "figures"
SEED = 20260828
QUANTILE = 0.10


def _canon(value):
    try:
        return Composition(str(value)).reduced_formula
    except Exception:
        return None


def build_unique_structure_cohort():
    meta = pd.read_parquet(PROC / "views_meta.parquet").copy()
    soap = np.load(PROC / "soap_geo.npy")
    if len(meta) != len(soap):
        raise ValueError("SOAP and metadata rows are not aligned")
    meta["canon"] = meta["formula"].map(_canon)
    counts = meta["canon"].value_counts(dropna=True)
    unique = meta["canon"].map(counts).eq(1)
    meta = meta.loc[unique].copy()
    meta["soap_row"] = meta.index.to_numpy()

    dual = pd.read_parquet(PROC / "dual_channel.parquet")
    needed = [
        "kappa_L_snyder_median", "kappa_e_n", "kappa_e_p", "PF_n", "PF_p",
        "T_electronic", "n_mp_polymorphs",
    ]
    valid = dual["n_mp_polymorphs"].eq(1)
    for col in needed[:-1]:
        valid &= np.isfinite(pd.to_numeric(dual[col], errors="coerce"))
    for col in ["kappa_L_snyder_median", "kappa_e_n", "kappa_e_p", "PF_n", "PF_p"]:
        valid &= pd.to_numeric(dual[col], errors="coerce") > 0
    dual = dual.loc[valid].copy()

    cohort = meta.merge(dual, on="canon", how="inner", validate="one_to_one")
    # Both tables carry lattice summary columns.  Keep the formula-level dual-channel
    # values under stable unsuffixed names for downstream diagnostics.
    for column in ["bulk_vrh", "shear_vrh", "debye", "density"]:
        dual_column = f"{column}_y"
        if dual_column in cohort:
            cohort[column] = cohort[dual_column]
    x = soap[cohort["soap_row"].to_numpy(int)]
    if not cohort["material_id"].is_unique or not cohort["canon"].is_unique:
        raise ValueError("Expected one structure and one formula per row")
    return cohort.reset_index(drop=True), x


def _stable_order(ids):
    return np.argsort([
        int.from_bytes(hashlib.blake2b(str(x).encode(), digest_size=8).digest(), "little")
        for x in ids
    ])


def embed_structure(soap):
    scaled = StandardScaler().fit_transform(soap)
    pca = PCA(n_components=min(50, soap.shape[1]), random_state=SEED)
    z = pca.fit_transform(scaled)
    return z[:, :2], z, float(pca.explained_variance_ratio_[:2].sum())


def assign_groups(cohort, carrier):
    log_ke = np.log10(cohort[f"kappa_e_{carrier}"].to_numpy(float))
    log_kl = np.log10(cohort["kappa_L_snyder_median"].to_numpy(float))
    log_ze = np.log10(
        cohort[f"PF_{carrier}"].to_numpy(float)
        * cohort["T_electronic"].to_numpy(float)
        / cohort[f"kappa_e_{carrier}"].to_numpy(float)
    )
    log_pf = np.log10(cohort[f"PF_{carrier}"].to_numpy(float))
    low_ke_cut = np.quantile(log_ke, QUANTILE)
    low_kl_cut = np.quantile(log_kl, QUANTILE)
    high_ze_cut = np.quantile(log_ze, 1 - QUANTILE)
    high_pf_cut = np.quantile(log_pf, 1 - QUANTILE)
    return {
        "low_ke": log_ke <= low_ke_cut,
        "low_kl": log_kl <= low_kl_cut,
        "high_electronic_quality": log_ze >= high_ze_cut,
        "high_pf": log_pf >= high_pf_cut,
        "low_ke_cut_log10": low_ke_cut,
        "low_kl_cut_log10": low_kl_cut,
        "high_ze_cut_log10": high_ze_cut,
        "high_pf_cut_log10": high_pf_cut,
        "log_ke": log_ke,
        "log_kl": log_kl,
        "log_electronic_quality": log_ze,
        "log_pf": log_pf,
    }


def _category(a, b, a_label):
    out = np.full(len(a), "other", dtype=object)
    out[a & ~b] = a_label
    out[~a & b] = "low lattice kappa"
    out[a & b] = "both"
    return out


def separation_metrics(z, label, name, carrier, n_boot=1000):
    label = np.asarray(label, dtype=bool)
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED),
    )
    probability = cross_val_predict(model, z, label.astype(int), cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(label, probability)
    rng = np.random.default_rng(SEED + sum(map(ord, f"{carrier}-{name}")))
    auc_boot = []
    for _ in range(n_boot):
        ix = rng.integers(0, len(label), len(label))
        if np.unique(label[ix]).size == 2:
            auc_boot.append(roc_auc_score(label[ix], probability[ix]))
    auc_lo, auc_hi = np.percentile(auc_boot, [2.5, 97.5])

    nn = NearestNeighbors(n_neighbors=16, metric="euclidean").fit(z)
    neighbours = nn.kneighbors(return_distance=False)[:, 1:]
    local_positive = label[neighbours].mean(axis=1)
    prevalence = label.mean()
    enrichment = local_positive[label].mean() / prevalence
    positive_ix = np.flatnonzero(label)
    enrich_boot = np.empty(n_boot)
    for b in range(n_boot):
        ix = rng.choice(positive_ix, len(positive_ix), replace=True)
        enrich_boot[b] = local_positive[ix].mean() / prevalence
    enrich_lo, enrich_hi = np.percentile(enrich_boot, [2.5, 97.5])
    return {
        "carrier": carrier,
        "label": name,
        "n": len(label),
        "n_positive": int(label.sum()),
        "prevalence": prevalence,
        "structure_auc": auc,
        "auc_ci_lo": auc_lo,
        "auc_ci_hi": auc_hi,
        "knn15_enrichment": enrichment,
        "knn_ci_lo": enrich_lo,
        "knn_ci_hi": enrich_hi,
        "feature_space": "composition-blind SOAP, first 50 PCs",
    }


def experimental_low_kl(cohort):
    exp = pd.read_csv(PROC / "experimental_formula_targets.csv")
    joined = cohort[["canon"]].merge(exp[["canon", "kappa_L"]], on="canon", how="left")
    available = joined["kappa_L"].notna() & (joined["kappa_L"] > 0)
    low = np.zeros(len(cohort), dtype=bool)
    if available.any():
        cutoff = joined.loc[available, "kappa_L"].quantile(.25)
        low[available] = joined.loc[available, "kappa_L"] <= cutoff
    else:
        cutoff = np.nan
    return low, int(available.sum()), float(cutoff)


def plot_map(cohort, xy, groups, exp_low, explained):
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.3), sharex=True, sharey=True,
                             constrained_layout=True)
    schemes = [
        ("literal", "low_ke", "low electronic kappa"),
        ("thermoelectric", "high_electronic_quality", "high PF·T/κe"),
    ]
    colours = {
        "other": "#c8c8c8",
        "low electronic kappa": "#247ba0",
        "high PF·T/κe": "#159f8c",
        "low lattice kappa": "#f28e2b",
        "both": "#7b2cbf",
    }
    order = ["other", "low electronic kappa", "high PF·T/κe", "low lattice kappa", "both"]
    xlim = np.percentile(xy[:, 0], [.5, 99.5])
    ylim = np.percentile(xy[:, 1], [.5, 99.5])
    xpad = .03 * np.diff(xlim)[0]
    ypad = .03 * np.diff(ylim)[0]
    for row, carrier in enumerate(("n", "p")):
        for col, (scheme, key, label) in enumerate(schemes):
            ax = axes[row, col]
            category = _category(groups[carrier][key], groups[carrier]["low_kl"], label)
            for cat in order:
                mask = category == cat
                if not mask.any():
                    continue
                ax.scatter(
                    xy[mask, 0], xy[mask, 1],
                    s=6 if cat == "other" else 18,
                    alpha=.20 if cat == "other" else .82,
                    color=colours[cat], edgecolor="none", rasterized=True,
                    label=f"{cat} (N={mask.sum()})",
                    zorder=1 if cat == "other" else 2,
                )
            if exp_low.any():
                ax.scatter(xy[exp_low, 0], xy[exp_low, 1], s=75, marker="*",
                           facecolors="none", edgecolors="black", linewidths=.9,
                           label=f"experimental low κL (N={exp_low.sum()})", zorder=4)
            ax.set_title(
                f"{carrier}-type | "
                + ("literal low-κ channels" if scheme == "literal" else "thermoelectric channel quality")
            )
            ax.set_xlabel("geometry-SOAP PC1")
            ax.set_ylabel("geometry-SOAP PC2")
            ax.set_xlim(xlim[0] - xpad, xlim[1] + xpad)
            ax.set_ylim(ylim[0] - ypad, ylim[1] + ypad)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, fontsize=8, frameon=False, loc="best")
    fig.suptitle(
        f"Electronic- and lattice-favoured formulas in one target-free structure space\n"
        f"bottom/top deciles; grey = other; central 99% shown; "
        f"PC1+PC2 explain {explained:.1%} of SOAP variance",
        fontsize=14,
    )
    fig.savefig(FIG / "structure_space_channel_clusters.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_two_clusters(xy, groups, explained, quantile=.10, filename="structure_space_two_clusters.png"):
    """Minimal figure requested for presentation: high PF and low kappa_L only."""
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.7), sharex=True, sharey=True,
                             constrained_layout=True)
    xlim = np.percentile(xy[:, 0], [.5, 99.5])
    ylim = np.percentile(xy[:, 1], [.5, 99.5])
    xpad = .03 * np.diff(xlim)[0]
    ypad = .03 * np.diff(ylim)[0]
    for ax, carrier in zip(axes, ("n", "p")):
        high_pf = groups[carrier]["log_pf"] >= np.quantile(
            groups[carrier]["log_pf"], 1 - quantile
        )
        low_kl = groups[carrier]["log_kl"] <= np.quantile(
            groups[carrier]["log_kl"], quantile
        )
        both = high_pf & low_kl
        pf_only = high_pf & ~low_kl
        kl_only = low_kl & ~high_pf
        other = ~(high_pf | low_kl)

        ax.scatter(xy[other, 0], xy[other, 1], s=5, alpha=.13, color="#a9a9a9",
                   edgecolor="none", rasterized=True, label=f"other (N={other.sum()})", zorder=1)
        ax.scatter(xy[pf_only, 0], xy[pf_only, 1], s=28, alpha=.92, color="#00a6a6",
                   edgecolor="white", linewidth=.20, rasterized=True,
                   label=f"high PF only (N={pf_only.sum()})", zorder=3)
        ax.scatter(xy[kl_only, 0], xy[kl_only, 1], s=28, alpha=.92, color="#ff8c1a",
                   edgecolor="white", linewidth=.20, rasterized=True,
                   label=f"low κL only (N={kl_only.sum()})", zorder=3)
        ax.scatter(xy[both, 0], xy[both, 1], s=48, alpha=1.0, color="#7b2cbf",
                   edgecolor="black", linewidth=.35, rasterized=True,
                   label=f"high PF + low κL (N={both.sum()})", zorder=5)
        ax.set_xlim(xlim[0] - xpad, xlim[1] + xpad)
        ax.set_ylim(ylim[0] - ypad, ylim[1] + ypad)
        ax.set_xlabel("geometry-SOAP PC1")
        ax.set_ylabel("geometry-SOAP PC2")
        ax.set_title(f"{carrier}-type")
        ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle(
        "High electronic transport and low lattice conductivity in one structure space\n"
        f"high PF = top {quantile:.0%}; low κL = bottom {quantile:.0%}; grey = all others; "
        f"central 99% shown (PC1+PC2: {explained:.1%} variance)",
        fontsize=14,
    )
    fig.savefig(FIG / filename, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_top5_diagnostics(cohort, xy, groups):
    """Export the extreme intersection and group medians without calling it validated zT."""
    experimental = set(pd.read_csv(PROC / "experimental_formula_targets.csv")["canon"])
    candidate_rows, summary_rows = [], []
    properties = ["bulk_vrh", "shear_vrh", "debye", "density", "Eg_opt", "fe"]
    for carrier in ("n", "p"):
        g = groups[carrier]
        high_pf = g["log_pf"] >= np.quantile(g["log_pf"], .95)
        low_kl = g["log_kl"] <= np.quantile(g["log_kl"], .05)
        both = high_pf & low_kl
        quality = (
            cohort[f"PF_{carrier}"].to_numpy(float)
            * cohort["T_electronic"].to_numpy(float)
            / cohort[f"kappa_e_{carrier}"].to_numpy(float)
        )
        quality_percentile = pd.Series(quality).rank(pct=True).to_numpy() * 100
        for i in np.flatnonzero(both):
            row = cohort.iloc[i]
            elements = [el.symbol for el in Composition(row["canon"]).elements]
            candidate_rows.append({
                "carrier": carrier, "canon": row["canon"],
                "material_id": row["material_id"], "PF": row[f"PF_{carrier}"],
                "kappa_e_CRTA_scale": row[f"kappa_e_{carrier}"],
                "seebeck": row[f"S_{carrier}"],
                "sigma_CRTA_scale": row[f"sigma_{carrier}"],
                "kappa_L_snyder": row["kappa_L_snyder_median"],
                "PF_T_over_kappa_e_percentile": quality_percentile[i],
                "Eg_opt": row["Eg_opt"], "bulk_vrh": row["bulk_vrh"],
                "shear_vrh": row["shear_vrh"], "debye": row["debye"],
                "density": row["density"], "formation_energy": row["fe"],
                "contains_halogen": any(el in {"F", "Cl", "Br", "I"} for el in elements),
                "experimental_kappa_L_match": row["canon"] in experimental,
                "soap_pc1": xy[i, 0], "soap_pc2": xy[i, 1],
                "validation_status": "unvalidated_formula_level_proxy_candidate",
            })
        for label, mask in [
            ("all", np.ones(len(cohort), dtype=bool)),
            ("high_PF_top5", high_pf), ("low_kL_bottom5", low_kl),
            ("intersection", both),
        ]:
            record = {"carrier": carrier, "group": label, "n": int(mask.sum())}
            for prop in properties:
                record[f"median_{prop}"] = float(cohort.loc[mask, prop].median())
            summary_rows.append(record)
    pd.DataFrame(candidate_rows).to_csv(PROC / "structure_top5_joint_candidates.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(PROC / "structure_top5_group_summary.csv", index=False)


def main():
    FIG.mkdir(exist_ok=True)
    cohort, soap = build_unique_structure_cohort()
    order = _stable_order(cohort["canon"])
    cohort = cohort.iloc[order].reset_index(drop=True)
    soap = soap[order]
    xy, z, explained = embed_structure(soap)
    groups = {carrier: assign_groups(cohort, carrier) for carrier in ("n", "p")}
    exp_low, n_exp, exp_cut = experimental_low_kl(cohort)

    membership_rows = []
    metric_rows = []
    pair_rows = []
    for carrier in ("n", "p"):
        g = groups[carrier]
        literal = _category(g["low_ke"], g["low_kl"], "low electronic kappa")
        thermo = _category(g["high_electronic_quality"], g["low_kl"], "high PF*T/kappa_e")
        for i, row in cohort.iterrows():
            membership_rows.append({
                "canon": row["canon"], "material_id": row["material_id"],
                "carrier": carrier, "soap_pc1": xy[i, 0], "soap_pc2": xy[i, 1],
                "log10_kappa_e": g["log_ke"][i],
                "log10_kappa_L_snyder": g["log_kl"][i],
                "log10_PF_T_over_kappa_e": g["log_electronic_quality"][i],
                "low_kappa_e_decile": g["low_ke"][i],
                "low_kappa_L_decile": g["low_kl"][i],
                "high_electronic_quality_decile": g["high_electronic_quality"][i],
                "high_pf_decile": g["high_pf"][i],
                "literal_category": literal[i], "thermoelectric_category": thermo[i],
                "experimental_low_kappa_L_quartile": exp_low[i],
            })
        for label, values in [
            ("low electronic kappa", g["low_ke"]),
            ("high PF*T/kappa_e", g["high_electronic_quality"]),
            ("high power factor", g["high_pf"]),
            ("low Snyder lattice kappa", g["low_kl"]),
        ]:
            metric_rows.append(separation_metrics(z, values, label, carrier))
        for definition, electronic_good in [
            ("literal low kappa_e", g["low_ke"]),
            ("high PF*T/kappa_e", g["high_electronic_quality"]),
            ("high power factor", g["high_pf"]),
        ]:
            lattice_good = g["low_kl"]
            observed = int((electronic_good & lattice_good).sum())
            expected = float(electronic_good.sum() * lattice_good.sum() / len(cohort))
            union = int((electronic_good | lattice_good).sum())
            pair_rows.append({
                "carrier": carrier, "definition": definition, "n": len(cohort),
                "n_electronic_good": int(electronic_good.sum()),
                "n_lattice_good": int(lattice_good.sum()),
                "n_both": observed, "random_expected_both": expected,
                "overlap_enrichment": observed / expected,
                "jaccard": observed / union,
            })

    membership = pd.DataFrame(membership_rows)
    metrics = pd.DataFrame(metric_rows)
    pairs = pd.DataFrame(pair_rows)
    membership.to_csv(PROC / "structure_channel_membership.csv", index=False)
    metrics.to_csv(PROC / "structure_channel_separation.csv", index=False)
    pairs.to_csv(PROC / "structure_channel_pair_overlap.csv", index=False)
    pd.DataFrame([
        {"quantile": QUANTILE, "n_unique_structure_formulas": len(cohort),
         "soap_pc1_pc2_explained_variance": explained,
         "n_experimental_kappa_L_matches": n_exp,
         "experimental_low_quartile_cutoff": exp_cut}
    ]).to_csv(PROC / "structure_channel_map_metadata.csv", index=False)
    plot_map(cohort, xy, groups, exp_low, explained)
    plot_two_clusters(xy, groups, explained)
    plot_two_clusters(
        xy, groups, explained, quantile=.05,
        filename="structure_space_two_clusters_top5.png",
    )
    write_top5_diagnostics(cohort, xy, groups)
    print(f"unique formula/structure cohort: {len(cohort)}")
    print(f"experimental κL matches: {n_exp}; low-quartile stars: {exp_low.sum()}")
    print(metrics.to_string(index=False))
    print("\nelectronic/lattice good-set overlap")
    print(pairs.to_string(index=False))
    print("saved figures/structure_space_channel_clusters.png")
    print("saved figures/structure_space_two_clusters.png")
    print("saved figures/structure_space_two_clusters_top5.png")
    print("saved processed/structure_top5_joint_candidates.csv")


if __name__ == "__main__":
    main()
