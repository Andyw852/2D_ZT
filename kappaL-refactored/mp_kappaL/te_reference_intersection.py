"""Find the intersection of structure-like and electronic-like TE neighbourhoods.

This is deliberately a descriptor-similarity analysis, not a zT or transport
ranking.  It uses the same complete-case JARVIS cohort and the same joint map as
``unified_chemical_space.py``.  Fourteen predeclared thermoelectric benchmark
formula families are used only as reference landmarks.

For every non-reference material we compute two independent distances:

* structure/chemistry: composition + composition-blind SOAP -> 30 PCs;
* electronic structure: Eg, electron mass, hole mass, dielectric response.

The distance in each block is the mean distance to the three nearest benchmark
formula centroids.  The top 5% most similar materials in each block define two
neighbourhoods; their material-level intersection is highlighted in purple.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import hypergeom
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import config


COORD_PATH = config.PROC_DIR / "unified_chemical_space_coordinates.csv"
META_PATH = config.PROC_DIR / "jarvis_structure_electronic_targets.parquet"
COMP_PATH = config.PROC_DIR / "jarvis_composition_magpie_style.npy"
SOAP_PATH = config.PROC_DIR / "jarvis_geometry_soap.npy"
BENCHMARK_PATH = config.PROC_DIR / "unified_chemical_space_te_benchmarks.csv"

MEMBERSHIP_OUT = config.PROC_DIR / "te_reference_dual_space_membership.csv"
SUMMARY_OUT = config.PROC_DIR / "te_reference_dual_space_summary.json"
FIG_OUT = config.FIG_DIR / "te_reference_dual_space_intersection.png"

SEED = 20260828
STRUCTURE_PCS = 30
REFERENCE_NEIGHBOURS = 3
TOP_FRACTION = 0.05
ROBUSTNESS_CUTOFF = 0.80

COLORS = {
    "other": "#c6c8cc",
    "structure only": "#2f6fed",
    "electronic only": "#f28e2b",
    "intersection": "#8e44ad",
    "benchmark": "#00d8e8",
}


def mean_k_reference_distance(
    values: np.ndarray,
    reference_centroids: np.ndarray,
    k: int = REFERENCE_NEIGHBOURS,
) -> np.ndarray:
    """Mean Euclidean distance to the k nearest distinct reference families."""
    values = np.asarray(values, dtype=float)
    reference_centroids = np.asarray(reference_centroids, dtype=float)
    if values.ndim != 2 or reference_centroids.ndim != 2:
        raise ValueError("values and reference_centroids must be 2D")
    if values.shape[1] != reference_centroids.shape[1]:
        raise ValueError("values and reference centroids must share feature width")
    if not 1 <= k <= len(reference_centroids):
        raise ValueError("k must be between 1 and the number of reference families")
    distances = cdist(values, reference_centroids, metric="euclidean")
    nearest = np.partition(distances, kth=k - 1, axis=1)[:, :k]
    return nearest.mean(axis=1)


def empirical_similarity_percentile(distance: np.ndarray, pool: np.ndarray) -> np.ndarray:
    """Convert smaller-is-better distances to a 0--1 pool-relative percentile."""
    distance = np.asarray(distance, dtype=float)
    pool = np.asarray(pool, dtype=bool)
    reference = np.sort(distance[pool])
    if len(reference) == 0:
        raise ValueError("percentile reference pool is empty")
    return 1.0 - np.searchsorted(reference, distance, side="left") / len(reference)


def _formula_label(formula: str) -> str:
    return re.sub(r"(\d+)", r"$_{\1}$", str(formula))


def _load_blocks() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    coords = pd.read_csv(COORD_PATH)
    meta = pd.read_parquet(META_PATH).reset_index(names="source_row")
    meta = meta[
        (meta["gap_ev"] > .01)
        & (meta["m_electron"] > 0)
        & (meta["m_hole"] > 0)
        & (meta["epsilon_geo"] > 0)
    ]
    frame = coords[["row_id", "joint_axis_1", "joint_axis_2"]].merge(
        meta, on="row_id", how="inner", validate="one_to_one"
    )
    if len(frame) != len(coords):
        raise ValueError(f"Coordinate/metadata mismatch: {len(frame)} != {len(coords)}")

    composition = np.load(COMP_PATH)
    geometry = np.load(SOAP_PATH)
    rows = frame["source_row"].to_numpy(int)
    structure_raw = np.column_stack([composition[rows], geometry[rows]])
    structure_scaled = StandardScaler().fit_transform(structure_raw)
    structure = PCA(n_components=STRUCTURE_PCS, random_state=SEED).fit_transform(
        structure_scaled
    )
    electronic = StandardScaler().fit_transform(np.column_stack([
        np.log1p(frame["gap_ev"]),
        np.log10(frame["m_electron"]),
        np.log10(frame["m_hole"]),
        np.log10(frame["epsilon_geo"]),
    ]))
    return frame, structure, electronic


def _reference_centroids(
    frame: pd.DataFrame,
    structure: np.ndarray,
    electronic: np.ndarray,
) -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray]:
    benchmarks = pd.read_csv(BENCHMARK_PATH)
    formula_by_id = benchmarks.drop_duplicates("row_id").set_index("row_id")[
        "benchmark_formula"
    ]
    frame = frame.copy()
    frame["benchmark_formula"] = frame["row_id"].map(formula_by_id)
    names, structure_centroids, electronic_centroids = [], [], []
    for name, group in frame[frame["benchmark_formula"].notna()].groupby(
        "benchmark_formula", sort=True
    ):
        index = group.index.to_numpy(int)
        names.append(str(name))
        structure_centroids.append(structure[index].mean(axis=0))
        electronic_centroids.append(electronic[index].mean(axis=0))
    if len(names) < REFERENCE_NEIGHBOURS:
        raise ValueError("Too few benchmark formula families")
    return (
        frame,
        names,
        np.asarray(structure_centroids),
        np.asarray(electronic_centroids),
    )


def _membership_for_references(
    structure: np.ndarray,
    electronic: np.ndarray,
    structure_refs: np.ndarray,
    electronic_refs: np.ndarray,
    pool: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    structure_distance = mean_k_reference_distance(structure, structure_refs)
    electronic_distance = mean_k_reference_distance(electronic, electronic_refs)
    structure_cutoff = np.quantile(structure_distance[pool], TOP_FRACTION)
    electronic_cutoff = np.quantile(electronic_distance[pool], TOP_FRACTION)
    structure_like = pool & (structure_distance <= structure_cutoff)
    electronic_like = pool & (electronic_distance <= electronic_cutoff)
    return structure_distance, electronic_distance, structure_like, electronic_like


def _build_membership() -> tuple[pd.DataFrame, dict]:
    frame, structure, electronic = _load_blocks()
    frame, reference_names, structure_refs, electronic_refs = _reference_centroids(
        frame, structure, electronic
    )
    is_benchmark = frame["benchmark_formula"].notna().to_numpy()
    pool = ~is_benchmark
    structure_distance, electronic_distance, structure_like, electronic_like = (
        _membership_for_references(
            structure, electronic, structure_refs, electronic_refs, pool
        )
    )
    intersection = structure_like & electronic_like

    # Leave one benchmark formula family out.  This measures whether membership
    # is driven by a single anchor family rather than by the reference set.
    loo_count = np.zeros(len(frame), dtype=float)
    for held_out in range(len(reference_names)):
        keep = np.arange(len(reference_names)) != held_out
        _, _, s_like, e_like = _membership_for_references(
            structure, electronic, structure_refs[keep], electronic_refs[keep], pool
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
    category[is_benchmark] = "benchmark"

    out = frame[[
        "row_id", "formula", "chemical_system", "joint_axis_1", "joint_axis_2",
        "gap_ev", "m_electron", "m_hole", "epsilon_geo", "benchmark_formula",
    ]].copy()
    out["structure_reference_distance"] = structure_distance
    out["electronic_reference_distance"] = electronic_distance
    out["structure_similarity_percentile"] = structure_percentile
    out["electronic_similarity_percentile"] = electronic_percentile
    out["joint_similarity"] = joint_similarity
    out["structure_like_top5"] = structure_like
    out["electronic_like_top5"] = electronic_like
    out["intersection_top5"] = intersection
    out["intersection_loo_stability"] = loo_stability
    out["robust_intersection"] = intersection & (loo_stability >= ROBUSTNESS_CUTOFF)
    out["category"] = category

    n_pool = int(pool.sum())
    n_structure = int(structure_like.sum())
    n_electronic = int(electronic_like.sum())
    n_both = int(intersection.sum())
    expected = n_structure * n_electronic / n_pool
    p_value = float(hypergeom.sf(n_both - 1, n_pool, n_structure, n_electronic))
    summary = {
        "n_complete_case": len(out),
        "n_benchmark_entries": int(is_benchmark.sum()),
        "n_benchmark_formula_families": len(reference_names),
        "n_nonbenchmark_pool": n_pool,
        "top_fraction_each_space": TOP_FRACTION,
        "reference_neighbours": REFERENCE_NEIGHBOURS,
        "n_structure_like": n_structure,
        "n_electronic_like": n_electronic,
        "n_intersection": n_both,
        "n_robust_intersection_loo_ge_0_8": int(out["robust_intersection"].sum()),
        "random_expected_intersection": expected,
        "overlap_enrichment": n_both / expected,
        "hypergeometric_p_value": p_value,
        "reference_formulas": reference_names,
        "interpretation": (
            "descriptor similarity to predeclared TE references; not a zT, PF, "
            "weighted-mobility, or independent kappa_L ranking"
        ),
    }
    return out, summary


def _annotate_overlap(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    labels = overlap.sort_values(
        ["intersection_loo_stability", "joint_similarity"], ascending=False
    ).drop_duplicates("formula").head(6)
    labels = labels.sort_values("structure_similarity_percentile", ascending=False)
    text_y = np.linspace(.997, .957, len(labels))
    for y, row in zip(text_y, labels.itertuples(index=False)):
        ax.annotate(
            _formula_label(row.formula),
            xy=(row.electronic_similarity_percentile, row.structure_similarity_percentile),
            xytext=(.906, y),
            textcoords="data",
            fontsize=8.2,
            ha="left",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": .55},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .84, "pad": .45},
            zorder=8,
        )


def _plot(frame: pd.DataFrame, summary: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.3, 6.7), constrained_layout=True)
    order = ["other", "structure only", "electronic only", "intersection"]
    sizes = {"other": 7, "structure only": 26, "electronic only": 26, "intersection": 52}
    alpha = {"other": .16, "structure only": .84, "electronic only": .84, "intersection": .96}
    labels = {
        "other": "other",
        "structure only": "structure-like only",
        "electronic only": "electronic-like only",
        "intersection": "intersection",
    }

    ax = axes[0]
    for category in order:
        subset = frame[frame["category"] == category]
        ax.scatter(
            subset["joint_axis_1"], subset["joint_axis_2"],
            s=sizes[category], color=COLORS[category], alpha=alpha[category],
            edgecolors="#222222" if category == "intersection" else "none",
            linewidths=.45 if category == "intersection" else 0,
            rasterized=category == "other", zorder=1 if category == "other" else 3,
        )
    benchmark = frame[frame["category"] == "benchmark"]
    ax.scatter(
        benchmark["joint_axis_1"], benchmark["joint_axis_2"],
        marker="*", s=105, color=COLORS["benchmark"], edgecolors="#111111",
        linewidths=.7, zorder=6,
    )
    robust = frame[frame["robust_intersection"]]
    ax.scatter(
        robust["joint_axis_1"], robust["joint_axis_2"], s=92,
        facecolors="none", edgecolors="#111111", linewidths=.9, zorder=7,
    )
    xlim = np.percentile(frame["joint_axis_1"], [.5, 99.5])
    ylim = np.percentile(frame["joint_axis_2"], [.5, 99.5])
    ax.set_xlim(xlim[0] - .04 * np.ptp(xlim), xlim[1] + .04 * np.ptp(xlim))
    ax.set_ylim(ylim[0] - .04 * np.ptp(ylim), ylim[1] + .04 * np.ptp(ylim))
    ax.set_xlabel("joint chemical axis 1")
    ax.set_ylabel("joint chemical axis 2")
    ax.set_title("Where the two neighbourhoods meet\nin the joint structure–electronic map")

    ax = axes[1]
    ax.add_patch(Rectangle(
        (.95, .95), .05, .05, facecolor=COLORS["intersection"],
        edgecolor="none", alpha=.10, zorder=0,
    ))
    for category in order:
        subset = frame[frame["category"] == category]
        ax.scatter(
            subset["electronic_similarity_percentile"],
            subset["structure_similarity_percentile"],
            s=sizes[category] + (4 if category != "other" else 0),
            color=COLORS[category], alpha=alpha[category],
            edgecolors="#222222" if category == "intersection" else "none",
            linewidths=.45 if category == "intersection" else 0,
            rasterized=category == "other", zorder=1 if category == "other" else 3,
        )
    ax.scatter(
        benchmark["electronic_similarity_percentile"],
        benchmark["structure_similarity_percentile"],
        marker="*", s=105, color=COLORS["benchmark"], edgecolors="#111111",
        linewidths=.7, zorder=6,
    )
    ax.scatter(
        robust["electronic_similarity_percentile"],
        robust["structure_similarity_percentile"],
        s=92, facecolors="none", edgecolors="#111111", linewidths=.9, zorder=7,
    )
    ax.axvline(.95, color=COLORS["electronic only"], lw=1.1, ls="--")
    ax.axhline(.95, color=COLORS["structure only"], lw=1.1, ls="--")
    ax.set_xlim(.70, 1.005)
    ax.set_ylim(.70, 1.005)
    ax.set_xlabel("electronic-structure similarity percentile")
    ax.set_ylabel("structure/chemistry similarity percentile")
    ax.set_title("The intersection is explicit in score space\n(top 30% region shown)")
    ax.text(
        .705, .712,
        f"intersection: {summary['n_intersection']}  |  random expectation: "
        f"{summary['random_expected_intersection']:.1f}  |  enrichment: "
        f"{summary['overlap_enrichment']:.2f}×",
        ha="left", va="bottom", fontsize=9.5,
    )
    _annotate_overlap(ax, frame[frame["intersection_top5"]])

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6,
               markerfacecolor=COLORS["other"], markeredgecolor="none",
               label=f"other (N={(frame.category == 'other').sum():,})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor=COLORS["structure only"], markeredgecolor="none",
               label=f"structure-like only (N={(frame.category == 'structure only').sum()})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor=COLORS["electronic only"], markeredgecolor="none",
               label=f"electronic-like only (N={(frame.category == 'electronic only').sum()})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=8,
               markerfacecolor=COLORS["intersection"], markeredgecolor="#222222",
               label=f"intersection (N={summary['n_intersection']})"),
        Line2D([0], [0], marker="*", linestyle="none", markersize=10,
               markerfacecolor=COLORS["benchmark"], markeredgecolor="#111111",
               label=f"TE references (N={summary['n_benchmark_entries']})"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=9,
               markerfacecolor="none", markeredgecolor="#111111",
               label=f"robust overlap, LOO ≥80% (N={summary['n_robust_intersection_loo_ge_0_8']})"),
    ]
    axes[0].legend(
        handles=legend_handles, loc="lower left", ncol=2, frameon=True,
        framealpha=.88, fontsize=8.2,
    )
    fig.suptitle(
        "Intersection of structure-like and electronic-structure-like TE neighbourhoods\n"
        "Top 5% nearest to 14 benchmark-family centroids in each full descriptor block; "
        "benchmarks excluded from thresholds",
        fontsize=14,
    )
    fig.savefig(FIG_OUT, dpi=230, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for path in (COORD_PATH, META_PATH, COMP_PATH, SOAP_PATH, BENCHMARK_PATH):
        if not path.exists():
            raise FileNotFoundError(path)
    frame, summary = _build_membership()
    frame.to_csv(MEMBERSHIP_OUT, index=False)
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    _plot(frame, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nintersection formulas")
    print(
        frame[frame["intersection_top5"]][
            ["formula", "row_id", "joint_similarity", "intersection_loo_stability"]
        ].sort_values(
            ["intersection_loo_stability", "joint_similarity"], ascending=False
        ).to_string(index=False)
    )
    print(f"saved {MEMBERSHIP_OUT}")
    print(f"saved {SUMMARY_OUT}")
    print(f"saved {FIG_OUT}")


if __name__ == "__main__":
    main()
