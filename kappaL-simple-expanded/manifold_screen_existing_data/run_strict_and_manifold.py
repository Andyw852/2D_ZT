"""Replace weighted view averaging with an interpretable strict-AND rank manifold.

For every material pair, compute its within-view neighbour-rank percentile in
the structure and electronic spaces.  The joint dissimilarity is the worse of
the two ranks:

    r_AND(i, j) = max(r_structure(i, j), r_electronic(i, j)).

Thus excellent similarity in one view can never compensate for poor similarity
in the other.  PF, kL and zT are excluded from graph construction.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import pairwise_distances
from umap import UMAP


ROOT = Path(__file__).resolve().parent
BASE_SCRIPT = ROOT / "run_joint_manifold_screen.py"
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"
FIGURE_OUT = FIGURE_DIR / "strict_and_joint_manifold.png"
PDF_OUT = FIGURE_DIR / "strict_and_joint_manifold.pdf"
POINTS_OUT = OUTPUT_DIR / "strict_and_manifold_points.csv"
CANDIDATES_OUT = OUTPUT_DIR / "strict_and_candidate_ranking.csv"
VALIDATION_OUT = OUTPUT_DIR / "strict_and_sensitivity.csv"
SUMMARY_OUT = OUTPUT_DIR / "strict_and_summary.json"

SEED = 20260829
K_VALUES = (15, 30, 50)
MAIN_K = 30
N_HIGHLIGHT = 30


def load_base_module():
    spec = importlib.util.spec_from_file_location("joint_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(BASE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rank_percentile_matrices(
    structure: np.ndarray, electronic: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(structure)
    rows = np.arange(n)[:, None]

    def ranks(values: np.ndarray) -> np.ndarray:
        distance = pairwise_distances(values, metric="euclidean", n_jobs=-1).astype(
            np.float32
        )
        order = np.argsort(distance, axis=1)
        rank = np.empty(order.shape, dtype=np.int32)
        rank[rows, order] = np.arange(n, dtype=np.int32)[None, :]
        return rank

    structure_rank = ranks(structure)
    electronic_rank = ranks(electronic)
    worst_rank = np.maximum(structure_rank, electronic_rank)
    return structure_rank, electronic_rank, worst_rank


def strict_and_affinity(worst_rank: np.ndarray, k: int) -> tuple[sparse.csr_matrix, np.ndarray]:
    n = worst_rank.shape[0]
    work = worst_rank.copy()
    np.fill_diagonal(work, n + 1)
    indices = np.argpartition(work, kth=k - 1, axis=1)[:, :k]
    rows = np.repeat(np.arange(n), k)
    cols = indices.ravel()
    selected_rank = work[np.arange(n)[:, None], indices].astype(float) / (n - 1)
    local_threshold = np.maximum(selected_rank.max(axis=1), 1.0 / (n - 1))
    weights = np.exp(
        -np.square(selected_rank.ravel() / np.repeat(local_threshold, k))
    )
    graph = sparse.coo_matrix((weights, (rows, cols)), shape=(n, n)).tocsr()
    graph = graph.maximum(graph.T)
    graph.setdiag(0.0)
    graph.eliminate_zeros()
    return graph, local_threshold


def two_dimensional_layout(coordinates: np.ndarray) -> np.ndarray:
    return UMAP(
        n_neighbors=MAIN_K,
        min_dist=0.18,
        metric="euclidean",
        random_state=SEED,
        n_jobs=1,
    ).fit_transform(coordinates)


def select_independent_dual_candidates(frame: pd.DataFrame) -> pd.Index:
    pool = frame["unknown_to_local_zt_table"] & ~frame["seed_formula"]
    return frame.loc[pool, "dual_score"].nlargest(min(N_HIGHLIGHT, int(pool.sum()))).index


def strict_candidate_scores(base, frame: pd.DataFrame, coordinates: np.ndarray) -> pd.DataFrame:
    out = frame.copy()
    analog, distance, nearest_formula = base.score_from_seed_distance(coordinates, out)
    out["strict_and_analog_percentile"] = analog
    out["strict_and_distance_to_seed"] = distance
    out["strict_and_nearest_seed_formula"] = nearest_formula
    out["strict_and_dual_score"] = np.sqrt(
        out["strict_and_analog_percentile"] * out["dual_score"].clip(0.0, 1.0)
    )
    pool = out["unknown_to_local_zt_table"] & ~out["seed_formula"]
    out["strict_and_screen_percentile"] = np.nan
    out.loc[pool, "strict_and_screen_percentile"] = out.loc[
        pool, "strict_and_dual_score"
    ].rank(pct=True, method="average")
    return out


def validation_summary(base, frame: pd.DataFrame, coordinates_by_k: dict[int, np.ndarray]) -> pd.DataFrame:
    rows = []
    for k, coordinates in coordinates_by_k.items():
        fold = base.leave_one_formula_out({f"strict AND k={k}": coordinates}, frame)
        rows.append(fold)
    strict = pd.concat(rows, ignore_index=True)
    baseline = pd.read_csv(OUTPUT_DIR / "manifold_formula_retrieval.csv")
    keep = baseline[baseline["view"].isin(["structure S1", "electronic E2", "joint S1+E2"])]
    return pd.concat([keep, strict], ignore_index=True)


def summarize_validation(validation: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for view, group in validation.groupby("view", sort=False):
        values = group["retrieval_percentile"].to_numpy(float)
        boot = np.median(
            rng.choice(values, size=(5000, len(values)), replace=True), axis=1
        )
        rows.append(
            {
                "view": view,
                "median": float(np.median(values)),
                "median_ci_low": float(np.quantile(boot, 0.025)),
                "median_ci_high": float(np.quantile(boot, 0.975)),
                "recall": int(group["recovered_top10pct"].sum()),
                "total": int(group["held_formula"].nunique()),
            }
        )
    return pd.DataFrame(rows).set_index("view")


def annotate(base, ax, frame: pd.DataFrame, x: str, y: str) -> None:
    known = frame[frame["seed_formula"]].sort_values("external_zt_max", ascending=False)
    candidates = frame[frame["independent_dual_candidate"]].sort_values(
        "dual_score", ascending=False
    )
    base.annotate_formulas(ax, known, x, y, 7, 7.7)
    base.annotate_formulas(ax, candidates, x, y, 6, 7.5)


def plot(
    base,
    frame: pd.DataFrame,
    validation: pd.DataFrame,
    structure_rank: np.ndarray,
    electronic_rank: np.ndarray,
    worst_rank: np.ndarray,
    local_threshold: np.ndarray,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 12.0), constrained_layout=True)
    known = frame[frame["seed_formula"]]
    highlighted = frame[frame["independent_dual_candidate"]]
    other = frame[~frame["seed_formula"] & ~frame["independent_dual_candidate"]]

    for ax, prefix, title in (
        (axes[0, 0], "weighted", r"Weighted-sum joint: $d^2=(d_S^2+d_E^2)/2$"),
        (axes[0, 1], "strict", r"Strict AND joint: $r_{AND}=\max(r_S,r_E)$"),
    ):
        ax.scatter(
            other[f"{prefix}_x"], other[f"{prefix}_y"], s=8,
            color=base.COLORS["other"], alpha=0.22, edgecolors="none", rasterized=True,
        )
        ax.scatter(
            highlighted[f"{prefix}_x"], highlighted[f"{prefix}_y"], s=48,
            color=base.COLORS["candidate"], edgecolors="#222222", linewidths=0.45,
            zorder=5,
        )
        ax.scatter(
            known[f"{prefix}_x"], known[f"{prefix}_y"], s=88, marker="*",
            color=base.COLORS["known"], edgecolors="#111111", linewidths=0.65,
            zorder=7,
        )
        ax.set_xlabel("2D manifold layout 1")
        ax.set_ylabel("2D manifold layout 2")
        ax.set_title(title)
    annotate(base, axes[0, 1], frame, "strict_x", "strict_y")
    axes[0, 0].legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=base.COLORS["other"], markeredgecolor="none", markersize=6, label="other complete cases"),
            Line2D([0], [0], marker="*", linestyle="none", markerfacecolor=base.COLORS["known"], markeredgecolor="#111111", markersize=11, label="reported max zT≥1 formula match"),
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=base.COLORS["candidate"], markeredgecolor="#222222", markersize=8, label="top-30 dual-score candidates (manifold-independent)"),
        ],
        loc="best", fontsize=8.0, framealpha=0.88,
    )

    ax = axes[1, 0]
    anchor_rows = frame.index[frame["canon"].astype(str) == "SnSe"].tolist()
    anchor = anchor_rows[0] if anchor_rows else int(np.flatnonzero(frame["seed_formula"])[0])
    n = len(frame)
    xs = structure_rank[anchor].astype(float) / (n - 1)
    ys = electronic_rank[anchor].astype(float) / (n - 1)
    threshold = float(local_threshold[anchor])
    accepted = (worst_rank[anchor].astype(float) / (n - 1)) <= threshold + 1e-12
    accepted[anchor] = False
    ax.scatter(xs[~accepted], ys[~accepted], s=7, color=base.COLORS["other"], alpha=0.22, edgecolors="none", rasterized=True)
    ax.scatter(xs[accepted], ys[accepted], s=36, color=base.COLORS["candidate"], edgecolors="#222222", linewidths=0.4)
    ax.axvline(threshold, color=base.COLORS["structure"], lw=1.1, ls="--")
    ax.axhline(threshold, color=base.COLORS["electronic"], lw=1.1, ls="--")
    upper = min(1.0, max(0.12, threshold * 2.8))
    ax.set_xlim(-0.003, upper)
    ax.set_ylim(-0.003, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"structure neighbour-rank percentile $r_S$")
    ax.set_ylabel(r"electronic neighbour-rank percentile $r_E$")
    ax.set_title(
        f"Why AND is interpretable — anchor {frame.loc[anchor, 'canon']}\n"
        r"accepted neighbours minimize the worse rank; no cross-view compensation"
    )
    ax.text(
        threshold * 0.96, threshold * 0.96, f"AND neighbourhood\nmax rank ≤ {threshold:.3f}",
        ha="right", va="top", fontsize=8.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.4},
    )

    ax = axes[1, 1]
    order = [
        "structure S1", "electronic E2", "joint S1+E2",
        "strict AND k=15", "strict AND k=30", "strict AND k=50",
    ]
    grouped = summarize_validation(validation).loc[order]
    colors = [
        base.COLORS["structure"], base.COLORS["electronic"], base.COLORS["joint"],
        "#5b2c83", "#7d3c98", "#a569bd",
    ]
    ypos = np.arange(len(order))[::-1]
    bars = ax.barh(ypos, grouped["median"], color=colors, alpha=0.9)
    ax.errorbar(
        grouped["median"], ypos,
        xerr=np.vstack(
            [
                grouped["median"] - grouped["median_ci_low"],
                grouped["median_ci_high"] - grouped["median"],
            ]
        ),
        fmt="none", ecolor="#222222", elinewidth=0.9, capsize=2.5, zorder=5,
    )
    ax.axvline(0.5, color="#666666", lw=1.0, ls="--")
    ax.set_yticks(ypos, ["Structure S1", "Electronic E2", "Weighted joint", "Strict AND k=15", "Strict AND k=30", "Strict AND k=50"])
    ax.set_xlim(0.0, 1.04)
    ax.set_xlabel("held-out high-zT formula retrieval percentile")
    ax.set_title("Independent leave-one-formula-out test\nfull 12D diffusion coordinates; whiskers are formula-bootstrap 95% CI")
    for bar, (_, row) in zip(bars, grouped.iterrows()):
        ax.text(
            0.018,
            bar.get_y() + bar.get_height() / 2,
            f"{row['median']:.2f}; top-10% {int(row['recall'])}/{int(row['total'])}",
            va="center", ha="left", fontsize=8.3, color="white",
        )

    fig.suptitle(
        "Weight-free structure–electronic manifold: worst-view neighbour rank controls every edge\n"
        "Purple candidates are selected without manifold proximity, so their map position is not built into the colour",
        fontsize=14,
    )
    fig.savefig(FIGURE_OUT, dpi=230, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    base = load_base_module()
    descriptor = base.load_descriptor_module()
    frame = base.prepare_data(descriptor)
    structure = base.robust_block(frame, base.STRUCTURE_FEATURES)
    electronic = base.robust_block(frame, base.ELECTRONIC_REFINED_FEATURES)

    structure_rank, electronic_rank, worst_rank = rank_percentile_matrices(
        structure, electronic
    )
    strict_coordinates: dict[int, np.ndarray] = {}
    thresholds: dict[int, np.ndarray] = {}
    edge_counts: dict[int, int] = {}
    graph_components: dict[int, int] = {}
    for k in K_VALUES:
        graph, local_threshold = strict_and_affinity(worst_rank, k)
        strict_coordinates[k] = base.diffusion_coordinates(graph)
        thresholds[k] = local_threshold
        edge_counts[k] = int(graph.nnz // 2)
        graph_components[k] = int(sparse.csgraph.connected_components(graph, directed=False)[0])

    weighted_values = np.hstack(
        [structure / math.sqrt(2.0), electronic / math.sqrt(2.0)]
    )
    weighted_graph = base.knn_affinity(weighted_values, MAIN_K)
    weighted_coordinates = base.diffusion_coordinates(weighted_graph)
    weighted_layout = two_dimensional_layout(weighted_coordinates)
    strict_layout = two_dimensional_layout(strict_coordinates[MAIN_K])

    validation = validation_summary(base, frame, strict_coordinates)
    frame = strict_candidate_scores(base, frame, strict_coordinates[MAIN_K])
    independent = select_independent_dual_candidates(frame)
    frame["independent_dual_candidate"] = False
    frame.loc[independent, "independent_dual_candidate"] = True
    frame[["weighted_x", "weighted_y"]] = weighted_layout
    frame[["strict_x", "strict_y"]] = strict_layout

    candidates = frame[
        frame["unknown_to_local_zt_table"] & ~frame["seed_formula"]
    ].sort_values("strict_and_dual_score", ascending=False)
    candidate_columns = [
        "row_id", "formula", "canon", "chemical_system", "preferred_carrier",
        "strict_and_nearest_seed_formula", "strict_and_analog_percentile",
        "dual_score", "strict_and_dual_score", "strict_and_screen_percentile",
    ]
    candidates[candidate_columns].head(100).to_csv(CANDIDATES_OUT, index=False)
    validation.to_csv(VALIDATION_OUT, index=False)
    point_columns = [
        "row_id", "formula", "canon", "seed_formula", "external_zt_max",
        "unknown_to_local_zt_table", "independent_dual_candidate",
        "strict_and_nearest_seed_formula", "strict_and_analog_percentile",
        "dual_score", "strict_and_dual_score", "strict_and_screen_percentile",
        "weighted_x", "weighted_y", "strict_x", "strict_y",
    ]
    frame[point_columns].to_csv(POINTS_OUT, index=False)

    metrics = summarize_validation(validation)
    summary = {
        "no_new_first_principles_or_transport_calculation": True,
        "manifold_is_target_free": True,
        "n_materials": int(len(frame)),
        "n_seed_formulas": int(frame.loc[frame["seed_formula"], "canon"].nunique()),
        "joint_definition": "r_AND(i,j)=max(r_structure(i,j),r_electronic(i,j))",
        "view_weight_coefficients": None,
        "main_k": MAIN_K,
        "k_sensitivity": list(K_VALUES),
        "edge_counts": {str(k): value for k, value in edge_counts.items()},
        "graph_connected_components": {str(k): value for k, value in graph_components.items()},
        "metrics": metrics.to_dict(orient="index"),
        "visual_candidate_definition": "top 30 pre-existing dual scores; no manifold proximity used",
        "warning": "formula-level high-zT seeds are not phase-, doping-, or microstructure-resolved",
    }
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    plot(
        base, frame, validation, structure_rank, electronic_rank, worst_rank,
        thresholds[MAIN_K],
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nStrict-AND candidates")
    print(candidates[candidate_columns].head(30).to_string(index=False))


if __name__ == "__main__":
    main()
