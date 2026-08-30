"""Build an expanded, simple-parameter structure/electronic TE reference map.

This is an independent follow-up to ``kappaL-refactored``.  It deliberately
uses a small set of interpretable scalar parameters, but enlarges the cohort by
using the transport-complete portion of the 55,723-entry JARVIS-DFT snapshot.

The two similarity blocks are:

* structure/chemistry (13 scalars): elemental statistics, density, volume per
  atom, cell size and lattice-shape measures;
* electronic/transport (7 scalars): band gap, dielectric magnitude/anisotropy,
  n/p Seebeck coefficient and n/p conductivity at the JARVIS fixed condition
  (600 K, 1e20 cm^-3, constant relaxation time for conductivity).

The map is a reference-similarity visualization, not an experimental zT
prediction.  Benchmarks are excluded from the top-5% thresholds.
"""
from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from pymatgen.core import Composition, Element
from scipy.spatial.distance import cdist
from scipy.stats import hypergeom
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler, StandardScaler

warnings.filterwarnings(
    "ignore",
    message="No Pauling electronegativity for .*",
    category=UserWarning,
)


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
RAW_PATH = (
    WORKSPACE
    / "jarvis_2d_te_atlas"
    / "data"
    / "raw"
    / "external"
    / "jarvis_kl"
    / "jdft_3d-8-18-2021.json"
)
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"
MEMBERSHIP_OUT = OUTPUT_DIR / "expanded_simple_space_membership.csv"
SUMMARY_OUT = OUTPUT_DIR / "expanded_simple_space_summary.json"
FEATURE_OUT = OUTPUT_DIR / "expanded_simple_space_features.csv"
FIG_OUT = FIGURE_DIR / "te_reference_dual_space_intersection_expanded.png"
FIG_PDF_OUT = FIGURE_DIR / "te_reference_dual_space_intersection_expanded.pdf"

SEED = 20260829
TOP_FRACTION = 0.05
REFERENCE_NEIGHBOURS = 3
ROBUSTNESS_CUTOFF = 0.80

TE_BENCHMARK_FORMULAS = (
    "Bi2Te3",
    "Bi2Se3",
    "PbTe",
    "PbSe",
    "PbS",
    "SnSe",
    "GeTe",
    "SnTe",
    "Cu2Se",
    "Mg2Si",
    "SiGe",
    "TiNiSn",
    "ZrNiSn",
    "SrTiO3",
)

STRUCTURE_FEATURES = (
    "n_elements",
    "composition_entropy",
    "mean_atomic_number",
    "std_atomic_number",
    "log_mean_atomic_mass",
    "mass_coefficient_variation",
    "mean_electronegativity",
    "std_electronegativity",
    "log_density",
    "log_volume_per_atom",
    "log_n_atoms",
    "log_cell_anisotropy",
    "angle_distortion",
)

ELECTRONIC_FEATURES = (
    "log1p_gap_ev",
    "log_dielectric_geo",
    "log_dielectric_anisotropy",
    "seebeck_n_uv_per_k",
    "seebeck_p_uv_per_k",
    "log1p_conductivity_n_raw",
    "log1p_conductivity_p_raw",
)

REQUIRED_RAW_FIELDS = (
    "optb88vdw_bandgap",
    "epsx",
    "epsy",
    "epsz",
    "n-Seebeck",
    "p-Seebeck",
    "n-powerfact",
    "p-powerfact",
    "ncond",
    "pcond",
)

COLORS = {
    "other": "#c6c8cc",
    "structure only": "#2f6fed",
    "electronic only": "#f28e2b",
    "intersection": "#8e44ad",
    "benchmark": "#00d8e8",
}


def _finite_number(value: object) -> bool:
    if value is None or (isinstance(value, str) and value.lower() in {"", "na", "nan", "none"}):
        return False
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _reduced_formula(formula: str) -> str:
    return Composition(str(formula)).reduced_formula


def _formula_label(formula: str) -> str:
    return re.sub(r"(\d+)", r"$_{\1}$", str(formula))


def _weighted_statistics(comp: Composition, field: str) -> tuple[float, float]:
    values: list[float] = []
    weights: list[float] = []
    for symbol, amount in comp.get_el_amt_dict().items():
        element = Element(symbol)
        value = getattr(element, field)
        if value is None:
            continue
        values.append(float(value))
        weights.append(float(amount))
    if not values:
        return float("nan"), float("nan")
    array = np.asarray(values, dtype=float)
    weight = np.asarray(weights, dtype=float)
    weight /= weight.sum()
    mean = float(np.sum(weight * array))
    std = float(np.sqrt(np.sum(weight * (array - mean) ** 2)))
    return mean, std


def _composition_features(formula: str) -> dict[str, float]:
    comp = Composition(str(formula))
    amounts = np.asarray(list(comp.get_el_amt_dict().values()), dtype=float)
    fractions = amounts / amounts.sum()
    mean_z, std_z = _weighted_statistics(comp, "Z")
    mean_mass, std_mass = _weighted_statistics(comp, "atomic_mass")
    mean_x, std_x = _weighted_statistics(comp, "X")
    return {
        "n_elements": float(len(amounts)),
        "composition_entropy": float(-np.sum(fractions * np.log(fractions))),
        "mean_atomic_number": mean_z,
        "std_atomic_number": std_z,
        "log_mean_atomic_mass": float(np.log1p(mean_mass)),
        "mass_coefficient_variation": float(std_mass / max(mean_mass, 1e-12)),
        "mean_electronegativity": mean_x,
        "std_electronegativity": std_x,
    }


def _structure_features(record: dict) -> dict[str, float]:
    atoms = record["atoms"]
    lattice = np.asarray(atoms["lattice_mat"], dtype=float)
    lengths = np.linalg.norm(lattice, axis=1)
    angles = np.asarray(atoms.get("angles", (90.0, 90.0, 90.0)), dtype=float)
    n_atoms = int(record.get("nat") or len(atoms["elements"]))
    volume = abs(float(np.linalg.det(lattice)))
    density = float(record["density"])
    return {
        "log_density": float(np.log1p(max(density, 0.0))),
        "log_volume_per_atom": float(np.log1p(volume / n_atoms)),
        "log_n_atoms": float(np.log1p(n_atoms)),
        "log_cell_anisotropy": float(np.log(max(lengths) / max(min(lengths), 1e-12))),
        "angle_distortion": float(np.mean(np.abs(angles - 90.0)) / 90.0),
    }


def _electronic_features(record: dict) -> dict[str, float]:
    gap = float(record["optb88vdw_bandgap"])
    eps = np.asarray([record["epsx"], record["epsy"], record["epsz"]], dtype=float)
    eps_geo = float(np.prod(eps) ** (1.0 / 3.0))
    eps_anisotropy = float(eps.max() / max(eps.min(), 1e-12))
    return {
        "gap_ev": gap,
        "dielectric_geo": eps_geo,
        "dielectric_anisotropy": eps_anisotropy,
        "seebeck_n_uv_per_k": float(record["n-Seebeck"]),
        "seebeck_p_uv_per_k": float(record["p-Seebeck"]),
        "power_factor_n_raw": float(record["n-powerfact"]),
        "power_factor_p_raw": float(record["p-powerfact"]),
        "conductivity_n_raw": float(record["ncond"]),
        "conductivity_p_raw": float(record["pcond"]),
        "log1p_gap_ev": float(np.log1p(gap)),
        "log_dielectric_geo": float(np.log(max(eps_geo, 1e-12))),
        "log_dielectric_anisotropy": float(np.log(max(eps_anisotropy, 1.0))),
        "log1p_conductivity_n_raw": float(np.log1p(max(float(record["ncond"]), 0.0))),
        "log1p_conductivity_p_raw": float(np.log1p(max(float(record["pcond"]), 0.0))),
    }


def load_expanded_cohort() -> pd.DataFrame:
    if not RAW_PATH.exists():
        raise FileNotFoundError(RAW_PATH)
    records = json.load(open(RAW_PATH, encoding="utf-8"))
    rows: list[dict] = []
    seen: set[str] = set()
    for record in records:
        jid = str(record.get("jid", ""))
        if not jid or jid in seen:
            continue
        seen.add(jid)
        if not all(_finite_number(record.get(field)) for field in REQUIRED_RAW_FIELDS):
            continue
        gap = float(record["optb88vdw_bandgap"])
        if gap <= 0.01:
            continue
        eps = [float(record[field]) for field in ("epsx", "epsy", "epsz")]
        if min(eps) <= 0:
            continue
        try:
            row = {
                "row_id": jid,
                "formula": str(record["formula"]),
                "reduced_formula": _reduced_formula(str(record["formula"])),
                "chemical_system": "-".join(sorted(Composition(record["formula"]).chemical_system.split("-"))),
            }
            row.update(_composition_features(record["formula"]))
            row.update(_structure_features(record))
            row.update(_electronic_features(record))
        except (ValueError, TypeError, KeyError, ZeroDivisionError):
            continue
        rows.append(row)
    frame = pd.DataFrame(rows)
    if len(frame) < 8000:
        raise RuntimeError(f"Expanded cohort unexpectedly small: {len(frame)}")
    if frame["row_id"].duplicated().any():
        raise ValueError("Duplicate JARVIS IDs remain after loading")
    return frame.reset_index(drop=True)


def _scale_block(frame: pd.DataFrame, fields: tuple[str, ...]) -> np.ndarray:
    raw = frame.loc[:, fields].to_numpy(dtype=float)
    imputed = SimpleImputer(strategy="median").fit_transform(raw)
    robust = RobustScaler(quantile_range=(10.0, 90.0)).fit_transform(imputed)
    robust = np.clip(robust, -6.0, 6.0)
    return StandardScaler().fit_transform(robust)


def _unit_inertia(block: np.ndarray) -> np.ndarray:
    centered = block - block.mean(axis=0, keepdims=True)
    norm = float(np.linalg.norm(centered, ord="fro"))
    if norm == 0:
        raise ValueError("Cannot normalize a zero-inertia block")
    return centered / norm


def _joint_coordinates(structure: np.ndarray, electronic: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
    fused = np.column_stack([_unit_inertia(structure), _unit_inertia(electronic)])
    coords = PCA(n_components=2, random_state=SEED).fit_transform(fused)
    if np.corrcoef(coords[:, 0], frame["log_mean_atomic_mass"])[0, 1] < 0:
        coords[:, 0] *= -1
    if np.corrcoef(coords[:, 1], frame["log1p_gap_ev"])[0, 1] < 0:
        coords[:, 1] *= -1
    return coords


def _reference_centroids(
    frame: pd.DataFrame,
    structure: np.ndarray,
    electronic: np.ndarray,
) -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray]:
    lookup = {_reduced_formula(formula): formula for formula in TE_BENCHMARK_FORMULAS}
    tagged = frame.copy()
    tagged["benchmark_formula"] = tagged["reduced_formula"].map(lookup)
    names: list[str] = []
    structure_centroids: list[np.ndarray] = []
    electronic_centroids: list[np.ndarray] = []
    for name, group in tagged[tagged["benchmark_formula"].notna()].groupby(
        "benchmark_formula", sort=True
    ):
        index = group.index.to_numpy(dtype=int)
        names.append(str(name))
        structure_centroids.append(structure[index].mean(axis=0))
        electronic_centroids.append(electronic[index].mean(axis=0))
    if len(names) < REFERENCE_NEIGHBOURS:
        raise ValueError(f"Too few benchmark families in expanded cohort: {len(names)}")
    return tagged, names, np.asarray(structure_centroids), np.asarray(electronic_centroids)


def mean_k_reference_distance(values: np.ndarray, centroids: np.ndarray, k: int) -> np.ndarray:
    distances = cdist(values, centroids, metric="euclidean")
    nearest = np.partition(distances, kth=k - 1, axis=1)[:, :k]
    return nearest.mean(axis=1)


def _membership_for_references(
    structure: np.ndarray,
    electronic: np.ndarray,
    structure_refs: np.ndarray,
    electronic_refs: np.ndarray,
    pool: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = min(REFERENCE_NEIGHBOURS, len(structure_refs))
    structure_distance = mean_k_reference_distance(structure, structure_refs, k)
    electronic_distance = mean_k_reference_distance(electronic, electronic_refs, k)
    structure_cutoff = float(np.quantile(structure_distance[pool], TOP_FRACTION))
    electronic_cutoff = float(np.quantile(electronic_distance[pool], TOP_FRACTION))
    structure_like = pool & (structure_distance <= structure_cutoff)
    electronic_like = pool & (electronic_distance <= electronic_cutoff)
    return structure_distance, electronic_distance, structure_like, electronic_like


def empirical_similarity_percentile(distance: np.ndarray, pool: np.ndarray) -> np.ndarray:
    reference = np.sort(np.asarray(distance, dtype=float)[pool])
    return 1.0 - np.searchsorted(reference, distance, side="left") / len(reference)


def build_membership(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    structure = _scale_block(frame, STRUCTURE_FEATURES)
    electronic = _scale_block(frame, ELECTRONIC_FEATURES)
    coords = _joint_coordinates(structure, electronic, frame)
    frame = frame.copy()
    frame["joint_axis_1"] = coords[:, 0]
    frame["joint_axis_2"] = coords[:, 1]
    frame, reference_names, structure_refs, electronic_refs = _reference_centroids(
        frame, structure, electronic
    )
    benchmark = frame["benchmark_formula"].notna().to_numpy()
    pool = ~benchmark
    structure_distance, electronic_distance, structure_like, electronic_like = (
        _membership_for_references(
            structure, electronic, structure_refs, electronic_refs, pool
        )
    )
    intersection = structure_like & electronic_like

    loo_count = np.zeros(len(frame), dtype=float)
    for held_out in range(len(reference_names)):
        keep = np.arange(len(reference_names)) != held_out
        _, _, s_like, e_like = _membership_for_references(
            structure,
            electronic,
            structure_refs[keep],
            electronic_refs[keep],
            pool,
        )
        loo_count += s_like & e_like
    loo_stability = loo_count / len(reference_names)

    structure_percentile = empirical_similarity_percentile(structure_distance, pool)
    electronic_percentile = empirical_similarity_percentile(electronic_distance, pool)
    joint_similarity = np.sqrt(structure_percentile * electronic_percentile)

    category = np.full(len(frame), "other", dtype=object)
    category[structure_like & ~electronic_like] = "structure only"
    category[~structure_like & electronic_like] = "electronic only"
    category[intersection] = "intersection"
    category[benchmark] = "benchmark"

    frame["structure_reference_distance"] = structure_distance
    frame["electronic_reference_distance"] = electronic_distance
    frame["structure_similarity_percentile"] = structure_percentile
    frame["electronic_similarity_percentile"] = electronic_percentile
    frame["joint_similarity"] = joint_similarity
    frame["structure_like_top5"] = structure_like
    frame["electronic_like_top5"] = electronic_like
    frame["intersection_top5"] = intersection
    frame["intersection_loo_stability"] = loo_stability
    frame["robust_intersection"] = intersection & (loo_stability >= ROBUSTNESS_CUTOFF)
    frame["category"] = category

    n_pool = int(pool.sum())
    n_structure = int(structure_like.sum())
    n_electronic = int(electronic_like.sum())
    n_intersection = int(intersection.sum())
    expected = n_structure * n_electronic / n_pool
    p_value = float(hypergeom.sf(n_intersection - 1, n_pool, n_structure, n_electronic))
    summary = {
        "source": str(RAW_PATH),
        "n_raw_jarvis_entries": 55723,
        "n_expanded_complete_case": int(len(frame)),
        "n_benchmark_entries": int(benchmark.sum()),
        "n_benchmark_formula_families": int(len(reference_names)),
        "n_nonbenchmark_pool": n_pool,
        "top_fraction_each_space": TOP_FRACTION,
        "reference_neighbours": REFERENCE_NEIGHBOURS,
        "n_structure_like": n_structure,
        "n_electronic_like": n_electronic,
        "n_intersection": n_intersection,
        "n_unique_intersection_formulas": int(frame.loc[intersection, "reduced_formula"].nunique()),
        "n_robust_intersection_loo_ge_0_8": int(frame["robust_intersection"].sum()),
        "random_expected_intersection": float(expected),
        "overlap_enrichment": float(n_intersection / expected),
        "hypergeometric_p_value": p_value,
        "reference_formulas_present": reference_names,
        "structure_features": list(STRUCTURE_FEATURES),
        "electronic_features": list(ELECTRONIC_FEATURES),
        "transport_condition": "JARVIS scalar transport: 600 K, 1e20 cm^-3, constant-tau conductivity",
        "interpretation": (
            "similarity to predeclared TE reference families in simple scalar "
            "structure/chemistry and electronic/transport blocks; not experimental zT"
        ),
    }
    return frame, summary


def _annotate_overlap(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    labels = (
        overlap.sort_values(
            ["intersection_loo_stability", "joint_similarity"], ascending=False
        )
        .drop_duplicates("reduced_formula")
        .head(8)
        .sort_values("structure_similarity_percentile", ascending=False)
    )
    text_y = np.linspace(0.998, 0.948, len(labels))
    for y, row in zip(text_y, labels.itertuples(index=False)):
        ax.annotate(
            _formula_label(row.reduced_formula),
            xy=(row.electronic_similarity_percentile, row.structure_similarity_percentile),
            xytext=(0.932, y),
            textcoords="data",
            fontsize=8.0,
            ha="right",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.55},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 0.45},
            zorder=8,
        )


def plot(frame: pd.DataFrame, summary: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.9), constrained_layout=True)
    order = ["other", "structure only", "electronic only", "intersection"]
    sizes = {"other": 5, "structure only": 20, "electronic only": 20, "intersection": 46}
    alpha = {"other": 0.13, "structure only": 0.80, "electronic only": 0.80, "intersection": 0.96}

    ax = axes[0]
    for category in order:
        subset = frame[frame["category"] == category]
        ax.scatter(
            subset["joint_axis_1"],
            subset["joint_axis_2"],
            s=sizes[category],
            color=COLORS[category],
            alpha=alpha[category],
            edgecolors="#222222" if category == "intersection" else "none",
            linewidths=0.4 if category == "intersection" else 0,
            rasterized=category == "other",
            zorder=1 if category == "other" else 3,
        )
    benchmark = frame[frame["category"] == "benchmark"]
    robust = frame[frame["robust_intersection"]]
    ax.scatter(
        benchmark["joint_axis_1"],
        benchmark["joint_axis_2"],
        marker="*",
        s=100,
        color=COLORS["benchmark"],
        edgecolors="#111111",
        linewidths=0.7,
        zorder=6,
    )
    ax.scatter(
        robust["joint_axis_1"],
        robust["joint_axis_2"],
        s=84,
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.9,
        zorder=7,
    )
    xlim = np.percentile(frame["joint_axis_1"], [0.25, 99.75])
    ylim = np.percentile(frame["joint_axis_2"], [0.25, 99.75])
    ax.set_xlim(xlim[0] - 0.04 * np.ptp(xlim), xlim[1] + 0.04 * np.ptp(xlim))
    ax.set_ylim(ylim[0] - 0.04 * np.ptp(ylim), ylim[1] + 0.04 * np.ptp(ylim))
    ax.set_xlabel("simple joint chemical axis 1")
    ax.set_ylabel("simple joint chemical axis 2")
    ax.set_title("Where the two neighbourhoods meet\nin the expanded simple-parameter map")

    ax = axes[1]
    ax.add_patch(
        Rectangle((0.95, 0.95), 0.05, 0.05, facecolor=COLORS["intersection"], edgecolor="none", alpha=0.10)
    )
    for category in order:
        subset = frame[frame["category"] == category]
        ax.scatter(
            subset["electronic_similarity_percentile"],
            subset["structure_similarity_percentile"],
            s=sizes[category] + (4 if category != "other" else 0),
            color=COLORS[category],
            alpha=alpha[category],
            edgecolors="#222222" if category == "intersection" else "none",
            linewidths=0.4 if category == "intersection" else 0,
            rasterized=category == "other",
            zorder=1 if category == "other" else 3,
        )
    ax.scatter(
        benchmark["electronic_similarity_percentile"],
        benchmark["structure_similarity_percentile"],
        marker="*",
        s=100,
        color=COLORS["benchmark"],
        edgecolors="#111111",
        linewidths=0.7,
        zorder=6,
    )
    ax.scatter(
        robust["electronic_similarity_percentile"],
        robust["structure_similarity_percentile"],
        s=84,
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.9,
        zorder=7,
    )
    ax.axvline(0.95, color=COLORS["electronic only"], lw=1.1, ls="--")
    ax.axhline(0.95, color=COLORS["structure only"], lw=1.1, ls="--")
    ax.set_xlim(0.70, 1.005)
    ax.set_ylim(0.70, 1.005)
    ax.set_xlabel("electronic/transport similarity percentile")
    ax.set_ylabel("structure/chemistry similarity percentile")
    ax.set_title("Explicit intersection in score space\n(top 30% region shown)")
    ax.text(
        0.705,
        0.712,
        f"intersection: {summary['n_intersection']}  |  random expectation: "
        f"{summary['random_expected_intersection']:.1f}  |  enrichment: "
        f"{summary['overlap_enrichment']:.2f}×",
        ha="left",
        va="bottom",
        fontsize=9.3,
    )
    _annotate_overlap(ax, frame[frame["intersection_top5"]])

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6, markerfacecolor=COLORS["other"], markeredgecolor="none", label=f"other (N={(frame.category == 'other').sum():,})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7, markerfacecolor=COLORS["structure only"], markeredgecolor="none", label=f"structure-like only (N={(frame.category == 'structure only').sum()})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7, markerfacecolor=COLORS["electronic only"], markeredgecolor="none", label=f"electronic/transport-like only (N={(frame.category == 'electronic only').sum()})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=8, markerfacecolor=COLORS["intersection"], markeredgecolor="#222222", label=f"intersection (N={summary['n_intersection']})"),
        Line2D([0], [0], marker="*", linestyle="none", markersize=10, markerfacecolor=COLORS["benchmark"], markeredgecolor="#111111", label=f"TE references (N={summary['n_benchmark_entries']})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=9, markerfacecolor="none", markeredgecolor="#111111", label=f"robust overlap, LOO ≥80% (N={summary['n_robust_intersection_loo_ge_0_8']})"),
    ]
    axes[0].legend(handles=legend_handles, loc="lower left", ncol=2, frameon=True, framealpha=0.88, fontsize=8.0)
    fig.suptitle(
        "Expanded intersection of simple structure and electronic/transport TE neighbourhoods\n"
        f"N={summary['n_expanded_complete_case']:,}; top 5% nearest to "
        f"{summary['n_benchmark_formula_families']} benchmark-family centroids in each block",
        fontsize=14,
    )
    fig.savefig(FIG_OUT, dpi=230, bbox_inches="tight")
    fig.savefig(FIG_PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    cohort = load_expanded_cohort()
    membership, summary = build_membership(cohort)

    feature_columns = [
        "row_id",
        "formula",
        "reduced_formula",
        "chemical_system",
        *STRUCTURE_FEATURES,
        "gap_ev",
        "dielectric_geo",
        "dielectric_anisotropy",
        "seebeck_n_uv_per_k",
        "seebeck_p_uv_per_k",
        "power_factor_n_raw",
        "power_factor_p_raw",
        "conductivity_n_raw",
        "conductivity_p_raw",
        *[field for field in ELECTRONIC_FEATURES if field not in {"seebeck_n_uv_per_k", "seebeck_p_uv_per_k"}],
    ]
    membership.loc[:, list(dict.fromkeys(feature_columns))].to_csv(FEATURE_OUT, index=False)
    output_columns = [
        "row_id",
        "formula",
        "reduced_formula",
        "chemical_system",
        "joint_axis_1",
        "joint_axis_2",
        "gap_ev",
        "seebeck_n_uv_per_k",
        "seebeck_p_uv_per_k",
        "power_factor_n_raw",
        "power_factor_p_raw",
        "benchmark_formula",
        "structure_reference_distance",
        "electronic_reference_distance",
        "structure_similarity_percentile",
        "electronic_similarity_percentile",
        "joint_similarity",
        "structure_like_top5",
        "electronic_like_top5",
        "intersection_top5",
        "intersection_loo_stability",
        "robust_intersection",
        "category",
    ]
    membership.loc[:, output_columns].to_csv(MEMBERSHIP_OUT, index=False)
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    plot(membership, summary)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nintersection formulas")
    print(
        membership.loc[
            membership["intersection_top5"],
            [
                "reduced_formula",
                "row_id",
                "joint_similarity",
                "intersection_loo_stability",
                "power_factor_n_raw",
                "power_factor_p_raw",
            ],
        ]
        .sort_values(["intersection_loo_stability", "joint_similarity"], ascending=False)
        .to_string(index=False)
    )
    print(f"saved {MEMBERSHIP_OUT}")
    print(f"saved {SUMMARY_OUT}")
    print(f"saved {FIG_OUT}")


if __name__ == "__main__":
    main()
