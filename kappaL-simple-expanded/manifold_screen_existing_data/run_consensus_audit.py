"""Audit the apparent overlap of high-zT seeds and dual-channel candidates.

This script does not perform new DFT, BTE, or phonon calculations.  It makes
the candidate definition explicit in two separate layers:

1. transport desirability from the existing cross-validated low-kL and PF
   scores; and
2. direct, non-compensatory similarity to the *same* high-zT seed in both the
   structure and electronic descriptor views.

The high-zT labels remain reduced-formula matches.  They are therefore drawn
as hollow stars and must not be interpreted as phase-resolved labels.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
BASE_SCRIPT = ROOT / "run_joint_manifold_screen.py"
STRICT_SCRIPT = ROOT / "run_strict_and_manifold.py"
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"

POINTS_IN = OUTPUT_DIR / "strict_and_manifold_points.csv"
OLD_POINTS_IN = OUTPUT_DIR / "joint_manifold_points.csv"
POINTS_OUT = OUTPUT_DIR / "consensus_audit_points.csv"
CANDIDATES_OUT = OUTPUT_DIR / "consensus_candidates.csv"
SUMMARY_OUT = OUTPUT_DIR / "consensus_audit_summary.json"
FIGURE_OUT = FIGURE_DIR / "consensus_structure_electronic_audit.png"
PDF_OUT = FIGURE_DIR / "consensus_structure_electronic_audit.pdf"

DUAL_QUANTILE = 0.90
VIEW_NEIGHBOUR_CUTOFF = 0.90


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_same_seed_consensus(
    frame: pd.DataFrame,
    structure_rank: np.ndarray,
    electronic_rank: np.ndarray,
) -> pd.DataFrame:
    """Find the seed minimizing the worse of the two directed view ranks.

    For a labelled row, all rows with the same reduced formula are excluded so
    its displayed position is leave-one-formula-out rather than a trivial
    self-match.  Candidate rows are compared with all seed rows.
    """
    out = frame.copy()
    n = len(out)
    all_seed = np.flatnonzero(out["seed_formula"].to_numpy(bool))
    formulas = out["canon"].astype(str).to_numpy()
    selected_seed = np.full(n, -1, dtype=int)
    selected_structure_rank = np.full(n, np.nan)
    selected_electronic_rank = np.full(n, np.nan)

    for row in range(n):
        eligible = all_seed
        if bool(out.iloc[row]["seed_formula"]):
            eligible = all_seed[formulas[all_seed] != formulas[row]]
        if len(eligible) == 0:
            continue
        worse = np.maximum(
            structure_rank[row, eligible], electronic_rank[row, eligible]
        )
        seed = int(eligible[int(np.argmin(worse))])
        selected_seed[row] = seed
        selected_structure_rank[row] = structure_rank[row, seed] / (n - 1)
        selected_electronic_rank[row] = electronic_rank[row, seed] / (n - 1)

    out["same_seed_structure_similarity"] = 1.0 - selected_structure_rank
    out["same_seed_electronic_similarity"] = 1.0 - selected_electronic_rank
    out["same_seed_and_similarity"] = np.minimum(
        out["same_seed_structure_similarity"],
        out["same_seed_electronic_similarity"],
    )
    out["same_seed_formula"] = [
        formulas[index] if index >= 0 else "" for index in selected_seed
    ]
    out["same_seed_row_id"] = [
        str(out.iloc[index]["row_id"]) if index >= 0 else ""
        for index in selected_seed
    ]
    return out


def add_candidate_flags(frame: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    out = frame.copy()
    pool = out["unknown_to_local_zt_table"] & ~out["seed_formula"]
    dual_threshold = float(out.loc[pool, "dual_score"].quantile(DUAL_QUANTILE))
    out["high_dual_score"] = pool & (out["dual_score"] >= dual_threshold)
    out["strict_same_seed_neighbour"] = pool & (
        (out["same_seed_structure_similarity"] >= VIEW_NEIGHBOUR_CUTOFF)
        & (out["same_seed_electronic_similarity"] >= VIEW_NEIGHBOUR_CUTOFF)
    )
    out["consensus_candidate"] = (
        out["high_dual_score"] & out["strict_same_seed_neighbour"]
    )
    return out, dual_threshold


def plot(frame: pd.DataFrame, summary: dict) -> None:
    other = frame[~frame["seed_formula"] & ~frame["consensus_candidate"]]
    seeds = frame[frame["seed_formula"]]
    candidates = frame[frame["consensus_candidate"]].sort_values(
        ["same_seed_and_similarity", "dual_score"], ascending=False
    )

    grey = "#c8cbd0"
    cyan = "#00d8e8"
    purple = "#8e44ad"
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 12.0), constrained_layout=True)

    ax = axes[0, 0]
    ax.scatter(
        other["same_seed_structure_similarity"],
        other["same_seed_electronic_similarity"],
        s=8, color=grey, alpha=0.22, edgecolors="none", rasterized=True,
    )
    ax.scatter(
        candidates["same_seed_structure_similarity"],
        candidates["same_seed_electronic_similarity"],
        s=54, color=purple, edgecolors="#222222", linewidths=0.5, zorder=5,
    )
    ax.scatter(
        seeds["same_seed_structure_similarity"],
        seeds["same_seed_electronic_similarity"],
        s=92, marker="*", facecolors="none", edgecolors="#00a9b7",
        linewidths=1.1, zorder=7,
    )
    ax.axvline(VIEW_NEIGHBOUR_CUTOFF, color="#2f6fed", ls="--", lw=1.1)
    ax.axhline(VIEW_NEIGHBOUR_CUTOFF, color="#f28e2b", ls="--", lw=1.1)
    ax.fill_between(
        [VIEW_NEIGHBOUR_CUTOFF, 1.0], VIEW_NEIGHBOUR_CUTOFF, 1.0,
        color=purple, alpha=0.08,
    )
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(0.0, 1.01)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("structure similarity to the selected high-zT seed\n(1 - within-view neighbour-rank percentile)")
    ax.set_ylabel("electronic similarity to the same seed\n(1 - within-view neighbour-rank percentile)")
    ax.set_title("Direct same-seed consensus plane\nupper-right means close in both views; no view compensation")

    ax = axes[0, 1]
    ax.scatter(
        other["strict_x"], other["strict_y"], s=8, color=grey,
        alpha=0.22, edgecolors="none", rasterized=True,
    )
    for row in candidates.itertuples(index=False):
        seed = frame.loc[frame["row_id"].astype(str) == str(row.same_seed_row_id)]
        if len(seed) == 1:
            ax.plot(
                [row.strict_x, float(seed.iloc[0]["strict_x"])],
                [row.strict_y, float(seed.iloc[0]["strict_y"])],
                color=purple, alpha=0.30, lw=0.8, zorder=2,
            )
    ax.scatter(
        candidates["strict_x"], candidates["strict_y"], s=54,
        color=purple, edgecolors="#222222", linewidths=0.5, zorder=5,
    )
    ax.scatter(
        seeds["strict_x"], seeds["strict_y"], s=92, marker="*",
        facecolors="none", edgecolors="#00a9b7", linewidths=1.1, zorder=7,
    )
    ax.set_xlabel("2D strict-AND manifold layout 1")
    ax.set_ylabel("2D strict-AND manifold layout 2")
    ax.set_title("The same points on the strict-AND UMAP\nlines link candidates to their same-seed analogue")

    ax = axes[1, 0]
    ax.scatter(
        other["electronic_score_percentile"],
        other["structure_score_percentile"],
        s=8, color=grey, alpha=0.22, edgecolors="none", rasterized=True,
    )
    ax.scatter(
        candidates["electronic_score_percentile"],
        candidates["structure_score_percentile"],
        s=54, color=purple, edgecolors="#222222", linewidths=0.5, zorder=5,
    )
    ax.scatter(
        seeds["electronic_score_percentile"],
        seeds["structure_score_percentile"],
        s=92, marker="*", facecolors="none", edgecolors="#00a9b7",
        linewidths=1.1, zorder=7,
    )
    ax.axvline(summary["dual_score_top10_threshold_electronic_min"], color="#f28e2b", ls=":", lw=0.9)
    ax.axhline(summary["dual_score_top10_threshold_structure_min"], color="#2f6fed", ls=":", lw=0.9)
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(0.0, 1.01)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("existing electronic/PF score percentile (600 K, fixed carrier density)")
    ax.set_ylabel("existing predicted low-kL score percentile (near 300 K labels)")
    ax.set_title("Transport-score plane\npurple also passes the independent same-seed similarity test")

    ax = axes[1, 1]
    ax.axis("off")
    audit_text = (
        "Why the previous pictures disagree\n\n"
        f"Complete cases: {summary['n_materials']:,}\n"
        f"High-zT reduced formulas: {summary['n_seed_formulas']}  "
        f"(drawn on {summary['n_seed_rows']} JARVIS structures)\n"
        f"Top-30 old-vs-independent candidate overlap: {summary['old_independent_overlap']}/30\n"
        f"Independent top-30 dual candidates in direct AND top 10%: "
        f"{summary['independent_top30_in_direct_and_top10']}/30\n"
        f"Strict consensus candidates shown in purple: {summary['n_consensus_candidates']}\n\n"
        "Old purple definition\n"
        "  joint-manifold proximity was multiplied into the colour score;\n"
        "  visual proximity was therefore partly selected in advance.\n\n"
        "Current hollow cyan stars\n"
        "  are reduced-formula matches, not phase-resolved samples.\n"
        "  Examples include C: g=0.45; SnS2: 3L;\n"
        "  Si: nano-bulk Si(model); InP/GaP: Reference.\n\n"
        "Interpretation\n"
        "  Purple means: predicted desirable transport AND close to the same\n"
        "  high-zT formula in both descriptor views.  It is a screening rule,\n"
        "  not independent evidence of high zT."
    )
    ax.text(
        0.02, 0.98, audit_text, transform=ax.transAxes, va="top", ha="left",
        fontsize=10.2, linespacing=1.35,
        bbox={"facecolor": "#f7f7f7", "edgecolor": "#dddddd", "pad": 0.8},
    )

    axes[0, 0].legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=grey, markeredgecolor="none", markersize=6, label="other complete cases"),
            Line2D([0], [0], marker="*", linestyle="none", markerfacecolor="none", markeredgecolor="#00a9b7", markersize=11, label="high-zT reduced-formula match (phase ambiguous)"),
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=purple, markeredgecolor="#222222", markersize=8, label="strict consensus candidate"),
        ],
        loc="lower left", fontsize=8.0, framealpha=0.9,
    )
    fig.suptitle(
        "Structure-electronic consensus audit using existing data only\n"
        "Candidate desirability, same-seed similarity, and 2D visualization are kept as separate claims",
        fontsize=14,
    )
    fig.savefig(FIGURE_OUT, dpi=230, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    base = load_module("joint_base_for_audit", BASE_SCRIPT)
    strict = load_module("strict_base_for_audit", STRICT_SCRIPT)
    descriptor = base.load_descriptor_module()
    frame = base.prepare_data(descriptor)
    structure = base.robust_block(frame, base.STRUCTURE_FEATURES)
    electronic = base.robust_block(frame, base.ELECTRONIC_REFINED_FEATURES)
    structure_rank, electronic_rank, _ = strict.rank_percentile_matrices(
        structure, electronic
    )
    frame = direct_same_seed_consensus(frame, structure_rank, electronic_rank)

    layout = pd.read_csv(POINTS_IN)[["row_id", "strict_x", "strict_y"]]
    frame = frame.merge(layout, on="row_id", how="left", validate="one_to_one")
    if frame[["strict_x", "strict_y"]].isna().any().any():
        raise ValueError("Strict layout rows do not align with rebuilt complete cases")
    frame, dual_threshold = add_candidate_flags(frame)

    pool = frame["unknown_to_local_zt_table"] & ~frame["seed_formula"]
    independent_top30 = set(frame.loc[pool].nlargest(30, "dual_score")["row_id"])
    old = pd.read_csv(OLD_POINTS_IN)
    old_top30 = set(old.loc[old["top_manifold_candidate"].astype(bool), "row_id"])
    independent_rows = frame[frame["row_id"].isin(independent_top30)]
    consensus = frame[frame["consensus_candidate"]].copy()

    summary = {
        "no_new_first_principles_or_transport_calculation": True,
        "n_materials": int(len(frame)),
        "n_seed_formulas": int(frame.loc[frame["seed_formula"], "canon"].nunique()),
        "n_seed_rows": int(frame["seed_formula"].sum()),
        "seed_warning": "reduced-formula matches; phase, doping, dimensionality, microstructure and measurement conditions are unresolved",
        "dual_score_quantile": DUAL_QUANTILE,
        "dual_score_threshold": dual_threshold,
        "view_neighbour_similarity_cutoff": VIEW_NEIGHBOUR_CUTOFF,
        "same_seed_rule": "choose seed minimizing max(r_structure, r_electronic); require both similarities >= 0.90",
        "n_consensus_candidates": int(frame["consensus_candidate"].sum()),
        "old_independent_overlap": int(len(old_top30 & independent_top30)),
        "independent_top30_in_direct_and_top10": int(
            (independent_rows["same_seed_and_similarity"] >= VIEW_NEIGHBOUR_CUTOFF).sum()
        ),
        "dual_score_top10_threshold_structure_min": float(
            frame.loc[frame["high_dual_score"], "structure_score_percentile"].min()
        ),
        "dual_score_top10_threshold_electronic_min": float(
            frame.loc[frame["high_dual_score"], "electronic_score_percentile"].min()
        ),
    }

    point_columns = [
        "row_id", "formula", "canon", "seed_formula", "external_zt_max",
        "unknown_to_local_zt_table", "same_seed_formula", "same_seed_row_id",
        "same_seed_structure_similarity", "same_seed_electronic_similarity",
        "same_seed_and_similarity", "structure_score_percentile",
        "electronic_score_percentile", "dual_score", "high_dual_score",
        "strict_same_seed_neighbour", "consensus_candidate", "strict_x", "strict_y",
    ]
    frame[point_columns].to_csv(POINTS_OUT, index=False)
    consensus[point_columns].sort_values(
        ["same_seed_and_similarity", "dual_score"], ascending=False
    ).to_csv(CANDIDATES_OUT, index=False)
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    plot(frame, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nStrict consensus candidates")
    sorted_consensus = consensus.sort_values(
        ["same_seed_and_similarity", "dual_score"], ascending=False
    )
    print(
        sorted_consensus[
            [
                "row_id", "formula", "same_seed_formula", "dual_score",
                "same_seed_structure_similarity",
                "same_seed_electronic_similarity",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
