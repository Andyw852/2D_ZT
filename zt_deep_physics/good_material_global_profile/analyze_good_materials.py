"""Compare high-ZT experimental samples with the global sample population.

The analysis is intentionally pairwise-complete: every property uses all peak-ZT
samples for which that property is available at the same sample's peak-ZT
temperature.  It never silently restricts the study to a tiny complete-case
intersection across all properties.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/zt_good_profile_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE.parent / "empirical" / "outputs" / "experimental_ZT_with_structure_metadata.csv"
DEFAULT_OUTPUT = HERE


@dataclass(frozen=True)
class Feature:
    column: str
    label: str
    unit: str
    category: str
    log_scale: bool = False


FEATURES = [
    Feature("temperature_peak_K", "Peak-ZT temperature", "K", "condition"),
    Feature("carrier_concentration_cm3", "Carrier concentration", "cm^-3", "electronic", True),
    Feature("abs_seebeck_uV_K", "Absolute Seebeck", "uV/K", "electronic"),
    Feature("sigma", "Electrical conductivity", "S/m", "electronic", True),
    Feature("power_factor_used_mW_mK2", "Power factor", "mW/mK^2", "electronic", True),
    Feature("kappa_lattice", "Lattice thermal conductivity", "W/mK", "thermal", True),
    Feature("kappa_total", "Total thermal conductivity", "W/mK", "thermal", True),
    Feature("kappa_electronic", "Electronic thermal conductivity", "W/mK", "thermal", True),
    Feature("mobility_cm2_Vs", "Mobility", "cm^2/Vs", "electronic", True),
    Feature("relative_density_pct", "Relative density", "%", "structure"),
    Feature("porosity_fraction", "Porosity fraction", "fraction", "structure"),
    Feature("grain_size_um", "Grain size", "um", "structure", True),
]

THRESHOLDS = [
    ("top_25pct", "Top 25%", "quantile", 0.75),
    ("top_10pct", "Top 10%", "quantile", 0.90),
    ("top_5pct", "Top 5%", "quantile", 0.95),
    ("zt_ge_1", "ZT >= 1", "absolute", 1.0),
]

# Match the source builder's accepted physical domain.  The power-factor guard
# is especially important because PF may be filled from S^2 sigma after the
# source columns were cleaned separately.
VALID_BOUNDS = {
    "temperature_peak_K": (100.0, 1500.0),
    "carrier_concentration_cm3": (1e14, 1e24),
    "abs_seebeck_uV_K": (0.0, 5000.0),
    "sigma": (1e-3, 1e9),
    "power_factor_used_mW_mK2": (0.0, 200.0),
    "kappa_lattice": (1e-3, 200.0),
    "kappa_total": (1e-3, 300.0),
    "kappa_electronic": (0.0, 200.0),
    "mobility_cm2_Vs": (1e-4, 1e5),
    "relative_density_pct": (20.0, 100.0),
    "porosity_fraction": (0.0, 0.8),
    "grain_size_um": (0.003, 500.0),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"sample_id", "zt_peak"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")
    frame = frame.drop_duplicates("sample_id", keep="first").copy()
    frame["zt_peak"] = pd.to_numeric(frame["zt_peak"], errors="coerce")
    frame = frame.loc[frame.zt_peak.notna() & frame.zt_peak.between(1e-8, 5.0)].copy()
    for feature in FEATURES:
        if feature.column in frame:
            frame[feature.column] = pd.to_numeric(frame[feature.column], errors="coerce")
            lo, hi = VALID_BOUNDS[feature.column]
            frame.loc[~frame[feature.column].between(lo, hi), feature.column] = np.nan
    # The empirical builder already applies source-level physical filters.  These
    # extra guards protect log transforms and parsed structure metadata.
    for feature in FEATURES:
        if feature.log_scale and feature.column in frame:
            frame.loc[frame[feature.column] <= 0, feature.column] = np.nan
    return frame


def quantiles(values: pd.Series, points=(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)) -> dict[str, float]:
    clean = values.dropna().astype(float)
    return {f"p{int(point * 100):02d}": float(clean.quantile(point)) for point in points}


def empirical_percentile(sorted_values: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Map values to the mid-rank empirical CDF of a reference population."""
    left = np.searchsorted(sorted_values, values, side="left")
    right = np.searchsorted(sorted_values, values, side="right")
    return 100.0 * (left + right) / (2.0 * len(sorted_values))


def evidence_grade(good_count: int, good_coverage: float) -> str:
    if good_count >= 500 and good_coverage >= 0.20:
        return "strong"
    if good_count >= 100 and good_coverage >= 0.05:
        return "moderate"
    return "limited"


def effect_label(auc_high_greater: float) -> tuple[str, float]:
    effect = 2.0 * abs(auc_high_greater - 0.5)
    if effect < 0.10:
        return "weak/mixed", effect
    direction = "higher" if auc_high_greater > 0.5 else "lower"
    strength = "strong" if effect >= 0.30 else "moderate"
    return f"{strength}; {direction}", effect


def weighted_stratified_auc(
    subset: pd.DataFrame, feature: str, threshold: float, stratum: pd.Series
) -> tuple[float, int]:
    """Probability that high-ZT ranks above other samples within strata."""
    work = subset[["zt_peak", feature]].copy()
    work["_stratum"] = stratum.reindex(work.index)
    work = work.dropna(subset=["_stratum"])
    numerator = 0.0
    denominator = 0
    used = 0
    for _, group in work.groupby("_stratum", observed=True):
        good = group.loc[group.zt_peak >= threshold, feature]
        other = group.loc[group.zt_peak < threshold, feature]
        if len(good) < 3 or len(other) < 3:
            continue
        comparisons = len(good) * len(other)
        numerator += mannwhitneyu(good, other, alternative="two-sided").statistic
        denominator += comparisons
        used += 1
    return (float(numerator / denominator), used) if denominator else (np.nan, 0)


def numeric_comparison(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    is_good = frame.zt_peak >= threshold
    temperature_strata = pd.qcut(frame.temperature_peak_K, q=5, duplicates="drop")
    family_strata = frame.material_family if "material_family" in frame else pd.Series(index=frame.index, dtype=object)
    rows = []
    for feature in FEATURES:
        if feature.column not in frame:
            continue
        subset = frame[["zt_peak", feature.column]].dropna()
        good = subset.loc[subset.zt_peak >= threshold, feature.column]
        other = subset.loc[subset.zt_peak < threshold, feature.column]
        if len(good) < 3 or len(other) < 3:
            continue
        global_q = quantiles(subset[feature.column])
        good_q = quantiles(good)
        u = mannwhitneyu(good, other, alternative="two-sided").statistic
        auc = float(u / (len(good) * len(other)))
        ks = ks_2samp(good, other, alternative="two-sided", method="auto")
        direction, effect = effect_label(auc)
        if feature.column == "temperature_peak_K":
            temp_auc, n_temp_strata = np.nan, 0
        else:
            temp_auc, n_temp_strata = weighted_stratified_auc(subset, feature.column, threshold, temperature_strata)
        family_auc, n_family_strata = weighted_stratified_auc(subset, feature.column, threshold, family_strata)
        core_lo, core_hi = good_q["p10"], good_q["p90"]
        global_in_core = subset[feature.column].between(core_lo, core_hi, inclusive="both").mean()
        high_in_global_core = good.between(global_q["p10"], global_q["p90"], inclusive="both").mean()
        good_coverage = len(good) / max(int(is_good.sum()), 1)
        row = {
            "feature": feature.column,
            "label": feature.label,
            "unit": feature.unit,
            "category": feature.category,
            "log_scale": feature.log_scale,
            "global_count": len(subset),
            "good_count": len(good),
            "global_coverage": len(subset) / len(frame),
            "good_coverage": good_coverage,
            **{f"global_{key}": value for key, value in global_q.items()},
            **{f"good_{key}": value for key, value in good_q.items()},
            "auc_probability_good_greater": auc,
            "rank_effect_strength": effect,
            "direction_and_strength": direction,
            "temperature_stratified_auc": temp_auc,
            "temperature_stratified_effect_strength": 2.0 * abs(temp_auc - 0.5) if np.isfinite(temp_auc) else np.nan,
            "temperature_strata_used": n_temp_strata,
            "family_stratified_auc": family_auc,
            "family_stratified_effect_strength": 2.0 * abs(family_auc - 0.5) if np.isfinite(family_auc) else np.nan,
            "family_strata_used": n_family_strata,
            "ks_statistic": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
            "global_fraction_in_good_p10_p90": float(global_in_core),
            "good_fraction_in_global_p10_p90": float(high_in_global_core),
            "evidence_grade": evidence_grade(len(good), good_coverage),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def build_enrichment_bins(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    for feature in FEATURES:
        if feature.column not in frame:
            continue
        subset = frame[["zt_peak", feature.column]].dropna().copy()
        if len(subset) < 30:
            continue
        subset["is_good"] = subset.zt_peak >= threshold
        baseline = subset.is_good.mean()
        try:
            subset["bin"] = pd.qcut(subset[feature.column], q=10, duplicates="drop")
        except ValueError:
            continue
        for index, (interval, group) in enumerate(subset.groupby("bin", observed=True), start=1):
            rate = group.is_good.mean()
            rows.append(
                {
                    "feature": feature.column,
                    "label": feature.label,
                    "unit": feature.unit,
                    "decile": index,
                    "lower": float(group[feature.column].min()),
                    "upper": float(group[feature.column].max()),
                    "n_total": len(group),
                    "n_good": int(group.is_good.sum()),
                    "good_rate": float(rate),
                    "paired_baseline_good_rate": float(baseline),
                    "enrichment": float(rate / baseline) if baseline > 0 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def threshold_sensitivity(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, label, kind, value in THRESHOLDS:
        threshold = float(frame.zt_peak.quantile(value)) if kind == "quantile" else float(value)
        good_mask = frame.zt_peak >= threshold
        for feature in FEATURES:
            if feature.column not in frame:
                continue
            subset = frame[["zt_peak", feature.column]].dropna()
            good = subset.loc[subset.zt_peak >= threshold, feature.column]
            other = subset.loc[subset.zt_peak < threshold, feature.column]
            if len(good) < 3 or len(other) < 3:
                continue
            q = quantiles(good, points=(0.10, 0.50, 0.90))
            auc = mannwhitneyu(good, other, alternative="two-sided").statistic / (len(good) * len(other))
            rows.append(
                {
                    "threshold_key": key,
                    "threshold_label": label,
                    "zt_threshold": threshold,
                    "n_good_all": int(good_mask.sum()),
                    "feature": feature.column,
                    "label": feature.label,
                    "unit": feature.unit,
                    "good_count": len(good),
                    "good_p10": q["p10"],
                    "good_p50": q["p50"],
                    "good_p90": q["p90"],
                    "auc_probability_good_greater": float(auc),
                }
            )
    return pd.DataFrame(rows)


def categorical_enrichment(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    for column, label in [("form", "Sample form"), ("material_family", "Material family")]:
        if column not in frame:
            continue
        subset = frame[["zt_peak", column]].dropna().copy()
        subset[column] = subset[column].astype(str).str.strip().replace("", np.nan)
        subset = subset.dropna(subset=[column])
        subset["is_good"] = subset.zt_peak >= threshold
        baseline = subset.is_good.mean()
        total_good = max(int(subset.is_good.sum()), 1)
        for value, group in subset.groupby(column):
            n_good = int(group.is_good.sum())
            rate = group.is_good.mean()
            rows.append(
                {
                    "feature": column,
                    "label": label,
                    "value": value,
                    "n_total": len(group),
                    "n_good": n_good,
                    "good_rate": float(rate),
                    "paired_baseline_good_rate": float(baseline),
                    "enrichment": float(rate / baseline) if baseline > 0 else np.nan,
                    "share_of_good": n_good / total_good,
                }
            )
    return pd.DataFrame(rows).sort_values(["feature", "enrichment"], ascending=[True, False])


def load_mechanism_priors() -> pd.DataFrame:
    """Load the project's scenario-model ranges without treating them as data."""
    path = HERE.parent / "outputs" / "top_decile_design_ranges.csv"
    if not path.exists():
        return pd.DataFrame(columns=["design_variable", "p10", "median", "p90", "evidence_scope"])
    priors = pd.read_csv(path)
    priors["evidence_scope"] = (
        "top 10% of sampled model scenarios after doping optimization; not a global material distribution"
    )
    return priors


def screening_ranges(numeric: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    rows = []
    robust_keys = {"top_10pct", "top_5pct", "zt_ge_1"}
    for row in numeric.itertuples(index=False):
        stable = sensitivity.loc[
            (sensitivity.feature == row.feature) & sensitivity.threshold_key.isin(robust_keys)
        ]
        robust_lo = float(stable.good_p10.max()) if len(stable) else np.nan
        robust_hi = float(stable.good_p90.min()) if len(stable) else np.nan
        if np.isfinite(robust_lo) and np.isfinite(robust_hi) and robust_lo > robust_hi:
            robust_lo, robust_hi = np.nan, np.nan
        if row.rank_effect_strength < 0.10:
            use = "no standalone cutoff; weak or mixed separation"
        elif row.evidence_grade == "limited":
            use = "hypothesis only; coverage is limited"
        else:
            use = "soft screening window; not a hard universal limit"
        if "higher" in row.direction_and_strength:
            constraint_kind = "soft_floor_retaining_90pct_of_paired_good"
            constraint_value = row.good_p10
        elif "lower" in row.direction_and_strength:
            constraint_kind = "soft_ceiling_retaining_90pct_of_paired_good"
            constraint_value = row.good_p90
        else:
            constraint_kind = "none"
            constraint_value = np.nan
        rows.append(
            {
                "feature": row.feature,
                "label": row.label,
                "unit": row.unit,
                "category": row.category,
                "typical_lower_p10": row.good_p10,
                "typical_upper_p90": row.good_p90,
                "broad_lower_p05": row.good_p05,
                "broad_upper_p95": row.good_p95,
                "robust_core_lower": robust_lo,
                "robust_core_upper": robust_hi,
                "direction_and_strength": row.direction_and_strength,
                "one_sided_constraint_kind": constraint_kind,
                "one_sided_constraint_value": constraint_value,
                "rank_effect_strength": row.rank_effect_strength,
                "good_coverage": row.good_coverage,
                "evidence_grade": row.evidence_grade,
                "recommended_use": use,
            }
        )
    return pd.DataFrame(rows)


def joint_rule_performance(frame: pd.DataFrame, threshold: float, numeric: pd.DataFrame) -> pd.DataFrame:
    lookup = numeric.set_index("feature")
    specifications = [
        (
            "PF soft floor",
            [("power_factor_used_mW_mK2", ">=", lookup.loc["power_factor_used_mW_mK2", "good_p10"])],
        ),
        (
            "kappa_total soft ceiling",
            [("kappa_total", "<=", lookup.loc["kappa_total", "good_p90"])],
        ),
        (
            "kappaL soft ceiling",
            [("kappa_lattice", "<=", lookup.loc["kappa_lattice", "good_p90"])],
        ),
        (
            "PF floor + kappa_total ceiling",
            [
                ("power_factor_used_mW_mK2", ">=", lookup.loc["power_factor_used_mW_mK2", "good_p10"]),
                ("kappa_total", "<=", lookup.loc["kappa_total", "good_p90"]),
            ],
        ),
        (
            "S-sigma window + kappa_total ceiling",
            [
                ("abs_seebeck_uV_K", ">=", lookup.loc["abs_seebeck_uV_K", "good_p10"]),
                ("abs_seebeck_uV_K", "<=", lookup.loc["abs_seebeck_uV_K", "good_p90"]),
                ("sigma", ">=", lookup.loc["sigma", "good_p10"]),
                ("sigma", "<=", lookup.loc["sigma", "good_p90"]),
                ("kappa_total", "<=", lookup.loc["kappa_total", "good_p90"]),
            ],
        ),
    ]
    rows = []
    for name, conditions in specifications:
        columns = sorted({condition[0] for condition in conditions})
        subset = frame.dropna(subset=columns).copy()
        passed = pd.Series(True, index=subset.index)
        descriptions = []
        for column, operator, value in conditions:
            passed &= subset[column].ge(value) if operator == ">=" else subset[column].le(value)
            descriptions.append(f"{column} {operator} {value:.8g}")
        good = subset.zt_peak >= threshold
        baseline = good.mean()
        n_pass = int(passed.sum())
        n_pass_good = int((passed & good).sum())
        precision = n_pass_good / n_pass if n_pass else np.nan
        retention = n_pass_good / max(int(good.sum()), 1)
        rows.append(
            {
                "rule": name,
                "conditions": " AND ".join(descriptions),
                "evaluable_n": len(subset),
                "evaluable_good_n": int(good.sum()),
                "pass_n": n_pass,
                "pass_good_n": n_pass_good,
                "global_pass_rate": float(passed.mean()),
                "high_zt_retention": retention,
                "high_zt_rate_among_pass": precision,
                "paired_baseline_high_zt_rate": float(baseline),
                "enrichment": float(precision / baseline) if baseline > 0 and np.isfinite(precision) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_ranges(numeric: pd.DataFrame, frame: pd.DataFrame, path: Path) -> None:
    records = []
    for row in numeric.itertuples(index=False):
        values = np.sort(frame[row.feature].dropna().to_numpy(float))
        points = np.array([row.good_p10, row.good_p50, row.good_p90])
        pct = empirical_percentile(values, points)
        records.append((row.label, *pct, row.evidence_grade))
    records.reverse()
    fig, ax = plt.subplots(figsize=(9.2, 7.2))
    ax.axvspan(10, 90, color="#e8eaed", alpha=0.8, label="global P10-P90")
    colors = {"strong": "#006d77", "moderate": "#e29578", "limited": "#8d99ae"}
    for y, (label, lo, med, hi, grade) in enumerate(records):
        ax.plot([lo, hi], [y, y], lw=5, solid_capstyle="round", color=colors[grade], alpha=0.9)
        ax.scatter(med, y, s=42, color="white", edgecolor=colors[grade], linewidth=1.8, zorder=3)
    ax.axvline(50, color="black", lw=0.8, ls="--")
    ax.set_yticks(range(len(records)), [record[0] for record in records])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Position of high-ZT P10 / median / P90 in the global percentile scale")
    ax.set_title("Where high-ZT samples sit relative to the global distribution")
    handles = [
        plt.Line2D([0], [0], color=colors[key], lw=5, label=f"{key} evidence")
        for key in ("strong", "moderate", "limited")
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_enrichment(enrichment: pd.DataFrame, path: Path) -> None:
    matrix = enrichment.pivot(index="label", columns="decile", values="enrichment")
    order = [feature.label for feature in FEATURES if feature.label in matrix.index]
    matrix = matrix.reindex(order)
    raw = matrix.to_numpy(float)
    values = np.full_like(raw, np.nan)
    valid = raw > 0
    values[valid] = np.log2(raw[valid])
    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    image = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_xticks(range(len(matrix.columns)), [f"D{int(value)}" for value in matrix.columns])
    ax.set_xlabel("Global property decile (D1 = lowest, D10 = highest)")
    ax.set_title("High-ZT enrichment across global property deciles")
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("log2(enrichment over pair-complete baseline)")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = matrix.iloc[i, j]
            if np.isfinite(value):
                color = "white" if abs(values[i, j]) > 1.1 else "black"
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=7, color=color)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_effect_coverage(numeric: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    colors = numeric.category.map(
        {"condition": "#6c757d", "electronic": "#0077b6", "thermal": "#d62828", "structure": "#2a9d8f"}
    )
    sizes = 35 + 85 * np.sqrt(numeric.good_count / numeric.good_count.max())
    ax.scatter(numeric.rank_effect_strength, 100 * numeric.good_coverage, s=sizes, c=colors, alpha=0.8)
    offsets = {
        "Carrier concentration": (-4, 8),
        "Mobility": (4, 10),
        "Electronic thermal conductivity": (4, -14),
        "Grain size": (4, 9),
        "Relative density": (6, 10),
        "Porosity fraction": (6, -14),
    }
    for row in numeric.itertuples(index=False):
        offset = offsets.get(row.label, (4, 3))
        ax.annotate(
            row.label,
            (row.rank_effect_strength, 100 * row.good_coverage),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.5,
        )
    ax.axvline(0.10, color="grey", lw=0.8, ls="--")
    ax.axvline(0.30, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("Rank separation strength: 2 x |AUC - 0.5|")
    ax.set_ylabel("Coverage among all high-ZT samples (%)")
    ax.set_xlim(0.08, max(0.72, numeric.rank_effect_strength.max() + 0.08))
    ax.set_title("A useful screening feature needs both separation and coverage")
    legend = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=8, label=label)
        for label, color in [("condition", "#6c757d"), ("electronic", "#0077b6"), ("thermal", "#d62828"), ("structure", "#2a9d8f")]
    ]
    ax.legend(handles=legend, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_robustness(sensitivity: pd.DataFrame, frame: pd.DataFrame, path: Path) -> None:
    selected = [
        "carrier_concentration_cm3",
        "abs_seebeck_uV_K",
        "sigma",
        "power_factor_used_mW_mK2",
        "kappa_lattice",
        "kappa_total",
    ]
    labels = {feature.column: feature.label for feature in FEATURES}
    threshold_order = [item[0] for item in THRESHOLDS]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharey=True)
    for ax, column in zip(axes.flat, selected):
        global_values = np.sort(frame[column].dropna().to_numpy(float))
        table = sensitivity.loc[sensitivity.feature.eq(column)].set_index("threshold_key").reindex(threshold_order)
        y = np.arange(len(table))
        lo = empirical_percentile(global_values, table.good_p10.to_numpy(float))
        med = empirical_percentile(global_values, table.good_p50.to_numpy(float))
        hi = empirical_percentile(global_values, table.good_p90.to_numpy(float))
        ax.axvspan(10, 90, color="#eceff1")
        for i in range(len(table)):
            ax.plot([lo[i], hi[i]], [y[i], y[i]], color="#264653", lw=4, solid_capstyle="round")
            ax.scatter(med[i], y[i], color="#e76f51", s=22, zorder=3)
        ax.set_title(labels[column])
        ax.set_xlim(0, 100)
        ax.set_yticks(y, table.threshold_label)
        ax.invert_yaxis()
        ax.set_xlabel("global percentile")
    fig.suptitle("High-performance ranges remain interpretable across ZT definitions", y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fmt_value(value: float, feature: str) -> str:
    if not np.isfinite(value):
        return "NA"
    if feature == "carrier_concentration_cm3" or abs(value) >= 1e5:
        return f"{value:.2e}"
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.3f}"


def chinese_direction(text: str) -> str:
    mapping = {
        "strong; higher": "明显偏高",
        "moderate; higher": "中等偏高",
        "strong; lower": "明显偏低",
        "moderate; lower": "中等偏低",
        "weak/mixed": "弱或混合",
    }
    return mapping.get(text, text)


def build_report(
    frame: pd.DataFrame,
    threshold: float,
    numeric: pd.DataFrame,
    screening: pd.DataFrame,
    enrichment: pd.DataFrame,
    categorical: pd.DataFrame,
    joint_rules: pd.DataFrame,
    mechanism_priors: pd.DataFrame,
    input_path: Path,
) -> str:
    n_good = int((frame.zt_peak >= threshold).sum())
    zt_ge_1 = int((frame.zt_peak >= 1.0).sum())
    lines = [
        "# 高性能热电材料相对全局的必要特征",
        "",
        "## 结论先行",
        "",
        f"本数据集共有 **{len(frame):,}** 个具有有效峰值 ZT 的实验样品。主分析把全局前 10% 定义为高性能组："
        f"`ZT_peak >= {threshold:.3f}`，共 **{n_good:,}** 个样品；另有 **{zt_ge_1:,}** 个样品满足 `ZT >= 1`。",
        "",
        "数据支持的核心不是一个孤立阈值，而是同时满足：**中等到较高的 |S|、足够高的电导/功率因子，以及较低的晶格或总热导**。"
        "载流子浓度给出可用工作窗，但单独区分高 ZT 的能力有限；迁移率和结构元数据当前覆盖不足或分离较弱，不能设成普适硬门槛。",
        "",
        "下表的 `P10-P90` 覆盖 80% 高性能配对样品，适合作为软筛选盒；它不是数学意义上的必要条件。",
        "",
        "| 性质 | 高性能组 N / 覆盖率 | 全局中位数 | 高性能典型 P10-P90 | 高性能中位数 | 相对全局 | 证据 |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    key_order = [
        "temperature_peak_K",
        "carrier_concentration_cm3",
        "abs_seebeck_uV_K",
        "sigma",
        "power_factor_used_mW_mK2",
        "kappa_lattice",
        "kappa_total",
        "kappa_electronic",
        "mobility_cm2_Vs",
        "relative_density_pct",
        "porosity_fraction",
        "grain_size_um",
    ]
    lookup = numeric.set_index("feature")
    for feature in key_order:
        if feature not in lookup.index:
            continue
        row = lookup.loc[feature]
        unit = f" {row.unit}" if row.unit else ""
        interval = f"{fmt_value(row.good_p10, feature)}-{fmt_value(row.good_p90, feature)}{unit}"
        lines.append(
            f"| {row.label} | {int(row.good_count):,} / {100*row.good_coverage:.1f}% | "
            f"{fmt_value(row.global_p50, feature)}{unit} | {interval} | "
            f"{fmt_value(row.good_p50, feature)}{unit} | {chinese_direction(row.direction_and_strength)} "
            f"(effect={row.rank_effect_strength:.2f}) | {row.evidence_grade} |"
        )

    useful = screening.loc[
        screening.evidence_grade.ne("limited") & screening.rank_effect_strength.ge(0.10)
    ].sort_values("rank_effect_strength", ascending=False)
    weak = screening.loc[screening.rank_effect_strength.lt(0.10)]
    lines += [
        "",
        "## 哪些特征最接近经验上的必要条件",
        "",
    ]
    if len(useful):
        for row in useful.itertuples(index=False):
            lines.append(
                f"- **{row.label}**：高性能组典型范围 "
                f"`{fmt_value(row.typical_lower_p10, row.feature)}-{fmt_value(row.typical_upper_p90, row.feature)} {row.unit}`；"
                f"{chinese_direction(row.direction_and_strength)}，高性能样品覆盖率 {100*row.good_coverage:.1f}%。"
            )
    if len(weak):
        weak_names = "、".join(weak.label.tolist())
        lines.append(f"- **不宜单独设阈值**：{weak_names} 的单变量秩分离效应小于 0.10。")

    lines += [
        "",
        "### 可直接用于初筛的单侧软约束",
        "",
        "下列界限各自保留该性质有数据的高性能样品约 90%；它们比双侧典型窗口更接近“必要条件”的操作定义。",
        "",
        "| 性质 | 软约束 | 配对高性能覆盖率 | 温度分层后效应 |",
        "|---|---:|---:|---:|",
    ]
    for row in useful.itertuples(index=False):
        if row.category == "condition":
            continue
        numeric_row = lookup.loc[row.feature]
        operator = ">=" if "floor" in row.one_sided_constraint_kind else "<="
        adjusted = numeric_row.temperature_stratified_effect_strength
        adjusted_text = f"{adjusted:.2f}" if np.isfinite(adjusted) else "NA"
        lines.append(
            f"| {row.label} | `{operator} {fmt_value(row.one_sided_constraint_value, row.feature)} {row.unit}` | "
            f"{100*row.good_coverage:.1f}% | {adjusted_text} |"
        )

    lines += [
        "",
        "### 联合规则的保留率与富集",
        "",
        "联合规则只在所需性质同时有数据的样品中评估；`保留率`越高越接近必要，`通过后高性能率/富集`越高越适合筛选。",
        "",
        "| 规则 | 可评估 N | 全局通过率 | 高性能保留率 | 通过后高性能率 | 富集 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in joint_rules.itertuples(index=False):
        lines.append(
            f"| {row.rule} | {row.evaluable_n:,} | {100*row.global_pass_rate:.1f}% | "
            f"{100*row.high_zt_retention:.1f}% | {100*row.high_zt_rate_among_pass:.1f}% | {row.enrichment:.2f}x |"
        )

    best_bins = (
        enrichment.sort_values(["feature", "enrichment"], ascending=[True, False])
        .groupby("feature", as_index=False)
        .first()
        .set_index("feature")
    )
    lines += [
        "",
        "## 富集而不是因果",
        "",
        "每个性质按全局十分位分箱后，最高富集区如下。富集倍数以该性质的成对完整样本为基线，因此不会把缺失率差异误写成整体 10% 基线。",
        "",
        "| 性质 | 最富集全局分位箱 | 数值范围 | 高性能率 | 相对基线富集 |",
        "|---|---:|---:|---:|---:|",
    ]
    for feature in key_order:
        if feature not in best_bins.index or feature not in lookup.index:
            continue
        row = best_bins.loc[feature]
        label = lookup.loc[feature, "label"]
        unit = lookup.loc[feature, "unit"]
        lines.append(
            f"| {label} | D{int(row.decile)} | {fmt_value(row.lower, feature)}-{fmt_value(row.upper, feature)} {unit} | "
            f"{100*row.good_rate:.1f}% | {row.enrichment:.2f}x |"
        )

    lines += [
        "",
        "## 结构与类别变量",
        "",
        "相对密度、孔隙率和晶粒尺寸来自文本元数据解析，覆盖远低于电子输运数据；样品形态和材料家族的富集还会混合研究热点、温区和发表选择，故只用于提出假设。",
        "",
    ]
    cat_show = categorical.loc[categorical.n_total.ge(30)].sort_values("enrichment", ascending=False).head(8)
    if len(cat_show):
        lines += [
            "| 类别 | 取值 | N | 高性能占比 | 富集 |",
            "|---|---|---:|---:|---:|",
        ]
        for row in cat_show.itertuples(index=False):
            lines.append(f"| {row.label} | {row.value} | {row.n_total:,} | {100*row.good_rate:.1f}% | {row.enrichment:.2f}x |")

    prior_meta = {
        "dos_mass_me": ("DOS effective mass", "m_e"),
        "conductivity_mass_me": ("Conductivity effective mass", "m_e"),
        "valley_degeneracy": ("Valley degeneracy", "count"),
        "elastic_modulus_2d_N_per_m": ("2D elastic modulus", "N/m"),
        "deformation_potential_eV": ("Deformation potential", "eV"),
        "porosity": ("Porosity", "fraction"),
        "group_velocity_m_s": ("Phonon group velocity", "m/s"),
        "lifetime_ps": ("Effective phonon lifetime", "ps"),
        "kappa_lattice_W_mK": ("Model lattice thermal conductivity", "W/mK"),
    }
    prior_show = mechanism_priors.loc[mechanism_priors.design_variable.isin(prior_meta)].copy()
    lines += [
        "",
        "## 与项目深层物理模型的对照",
        "",
        "实验全局表没有完整的有效质量、谷简并、形变势、二维弹性常数、群速度和声子寿命。"
        "下表来自项目原有 400 个透明参数情景中、各自优化掺杂后 ZT 前 10% 的 `P10-P90`，只能作为计算先验，不能称为相对全局的经验必要范围。",
        "",
        "| 深层变量 | 模型情景 P10-P90 | 中位数 | 证据性质 |",
        "|---|---:|---:|---|",
    ]
    for row in prior_show.itertuples(index=False):
        label, unit = prior_meta[row.design_variable]
        lines.append(
            f"| {label} | {fmt_value(row.p10, row.design_variable)}-{fmt_value(row.p90, row.design_variable)} {unit} | "
            f"{fmt_value(row.median, row.design_variable)} {unit} | 情景模型先验 |"
        )
    lines += [
        "",
        "特别需要保留的模型—实验张力：模型高表现情景允许较宽的孔隙率窗口，但实验结构元数据中的高 ZT 配对样品主要集中在高相对密度、低孔隙率。"
        "这说明‘用孔隙降 kappaL’不是免费的收益；真实材料中的电导损失、连通性和样品制备必须同时验证。",
        "",
        "模型与实验对 kappaL 的有利区间有明显重叠，但实验窗口更宽。对缺失的深层电子/声子变量，应把模型范围用于安排 DFT/声子计算优先级，而非直接淘汰材料。",
    ]

    lines += [
        "",
        "## 如何使用这些范围",
        "",
        "1. 初筛时先用 `screening_ranges.csv` 的 `P10-P90` 作为软筛选窗口，不要把边界外样品直接判死刑；",
        "2. 优先联合判断 `PF` 与 `kappa_total`（或 `kappaL + kappae`），因为 ZT 本身是这些量的组合；",
        "3. 用 `P05-P95` 检查候选是否明显偏离已知高性能经验域；",
        "4. 二维材料还需重新核对有效厚度定义，体单位的 sigma、kappa 和载流子浓度会随厚度约定改变；",
        "5. 对最终候选补充带隙、有效质量、谷简并、形变势、声子稳定性和三阶力常数。当前实验表并不完整包含这些深层变量。",
        "",
        "## 证据边界",
        "",
        "- 这是跨论文、跨材料家族、跨温度和跨制备状态的观察性比较，只能给出经验约束与优先级，不能证明单变量因果；",
        "- 峰值 ZT 温度本身是测试条件，不是材料固有常数；高温富集也受到测量温区覆盖影响；",
        "- 数值表同时给出温度五分位分层和材料家族内的秩效应；它们用于检查方向稳健性，不等同于完整因果校正；",
        "- 每个性质使用各自的最大配对样本，表中的 N 不同；严禁把小样本结构范围与万级电子输运范围视为同等证据；",
        "- PF 优先采用直接曲线，缺失时以同温度 `S^2 sigma` 补充；相对密度和晶粒尺寸是文本解析值；",
        "- P10-P90 按定义仍允许 20% 已知高性能样品位于区间外，因此这里称为软范围而非硬必要条件。",
        "",
        "## 可复现性",
        "",
        f"输入：`{input_path}`",
        "",
        "运行 `analyze_good_materials.py` 可重建全部表格、图和本报告。",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    outputs = output_dir / "outputs"
    figures = output_dir / "figures"
    outputs.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    frame = load_data(args.input.resolve())
    threshold = float(frame.zt_peak.quantile(0.90))
    numeric = numeric_comparison(frame, threshold)
    enrichment = build_enrichment_bins(frame, threshold)
    sensitivity = threshold_sensitivity(frame)
    categorical = categorical_enrichment(frame, threshold)
    mechanism_priors = load_mechanism_priors()
    ranges = screening_ranges(numeric, sensitivity)
    joint_rules = joint_rule_performance(frame, threshold, numeric)

    numeric.to_csv(outputs / "global_vs_good_numeric.csv", index=False)
    enrichment.to_csv(outputs / "quantile_enrichment_bins.csv", index=False)
    sensitivity.to_csv(outputs / "threshold_sensitivity.csv", index=False)
    categorical.to_csv(outputs / "categorical_enrichment.csv", index=False)
    ranges.to_csv(outputs / "screening_ranges.csv", index=False)
    joint_rules.to_csv(outputs / "joint_rule_performance.csv", index=False)
    mechanism_priors.to_csv(outputs / "mechanism_model_priors.csv", index=False)

    manifest = pd.DataFrame(
        [
            {"item": "input_file", "value": str(args.input.resolve())},
            {"item": "n_unique_peak_zt_samples", "value": len(frame)},
            {"item": "primary_definition", "value": "global top 10% peak ZT"},
            {"item": "primary_zt_threshold", "value": threshold},
            {"item": "n_primary_good", "value": int((frame.zt_peak >= threshold).sum())},
            {"item": "n_zt_ge_1", "value": int((frame.zt_peak >= 1.0).sum())},
            {"item": "pairing_rule", "value": "same sample, property interpolated at peak-ZT temperature"},
            {"item": "typical_window", "value": "P10-P90 of primary high-ZT group"},
            {"item": "broad_window", "value": "P05-P95 of primary high-ZT group"},
        ]
    )
    manifest.to_csv(outputs / "analysis_manifest.csv", index=False)

    plot_ranges(numeric, frame, figures / "01_good_ranges_on_global_percentiles.png")
    plot_enrichment(enrichment, figures / "02_decile_enrichment_heatmap.png")
    plot_effect_coverage(numeric, figures / "03_effect_vs_coverage.png")
    plot_robustness(sensitivity, frame, figures / "04_threshold_robustness.png")

    report = build_report(
        frame,
        threshold,
        numeric,
        ranges,
        enrichment,
        categorical,
        joint_rules,
        mechanism_priors,
        args.input.resolve(),
    )
    (outputs / "report.md").write_text(report, encoding="utf-8")
    print(f"Wrote analysis to {output_dir}")
    print(f"Primary high-ZT threshold: {threshold:.6g}; N={int((frame.zt_peak >= threshold).sum()):,}")


if __name__ == "__main__":
    main()
