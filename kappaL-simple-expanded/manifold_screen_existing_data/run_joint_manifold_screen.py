"""Screen thermoelectric analogues on a structure-electronic diffusion manifold.

Only existing local descriptors and labels are used.  No DFT, BTE, phonon or
new transport calculation is performed.  The manifold itself is target-free:
PF, kL and zT are never used to construct graph edges or diffusion coordinates.
Known starrydata2 high-zT formulas are weak seeds used only after construction.
"""
from __future__ import annotations

import importlib.util
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import RobustScaler
from umap import UMAP


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
DESCRIPTOR_SCRIPT = PARENT / "descriptor_fit_no_new_dft" / "fit_descriptor_space.py"
DUAL_SCORE_PATH = (
    PARENT
    / "descriptor_fit_no_new_dft"
    / "outputs"
    / "cross_validated_descriptor_space.csv"
)

OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"
POINTS_OUT = OUTPUT_DIR / "joint_manifold_points.csv"
CANDIDATES_OUT = OUTPUT_DIR / "manifold_candidate_ranking.csv"
VALIDATION_OUT = OUTPUT_DIR / "manifold_formula_retrieval.csv"
SUMMARY_OUT = OUTPUT_DIR / "joint_manifold_summary.json"
FIGURE_OUT = FIGURE_DIR / "joint_structure_electronic_manifold_screen.png"
PDF_OUT = FIGURE_DIR / "joint_structure_electronic_manifold_screen.pdf"

SEED = 20260829
K_GRAPH = 30
N_DIFFUSION = 12
N_CANDIDATES_PLOT = 30

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
ELECTRONIC_BASE_FEATURES = (
    "log1p_gap_ev",
    "log_dielectric_geo",
    "log_dielectric_anisotropy",
)
ELECTRONIC_REFINED_FEATURES = ELECTRONIC_BASE_FEATURES + (
    "log1p_mbj_gap",
    "log1p_electron_mass",
    "log1p_hole_mass",
    "log_mass_ratio",
    "log1p_electron_mass_spectral_ratio",
    "log1p_hole_mass_spectral_ratio",
    "log1p_electron_mass_complexity_proxy",
    "log1p_hole_mass_complexity_proxy",
)

COLORS = {
    "other": "#c8cbd0",
    "candidate": "#8e44ad",
    "known": "#00d8e8",
    "structure": "#2f6fed",
    "electronic": "#f28e2b",
    "joint": "#6f42c1",
}


def load_descriptor_module():
    spec = importlib.util.spec_from_file_location("descriptor_fit", DESCRIPTOR_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(DESCRIPTOR_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def formula_label(formula: str) -> str:
    return re.sub(r"(\d+)", r"$_{\1}$", str(formula))


def prepare_data(module) -> pd.DataFrame:
    raw = json.load(open(module.RAW_PATH, encoding="utf-8"))
    raw_by_jid = {str(record["jid"]): record for record in raw}
    frame = module.build_expanded_frame(raw_by_jid)
    frame = module.add_external_zt_labels(frame)
    dual = pd.read_csv(DUAL_SCORE_PATH)[
        [
            "row_id",
            "predicted_low_kL_score_raw",
            "predicted_electronic_score_raw",
            "structure_score_percentile",
            "electronic_score_percentile",
            "dual_score",
        ]
    ]
    frame = frame.merge(dual, on="row_id", how="left", validate="one_to_one")
    required = list(STRUCTURE_FEATURES + ELECTRONIC_REFINED_FEATURES)
    finite = np.isfinite(frame[required].to_numpy(float)).all(axis=1)
    frame = frame.loc[finite].reset_index(drop=True)
    frame["seed_formula"] = frame["external_high_zt_formula"].astype(bool)
    frame["unknown_to_local_zt_table"] = frame["external_zt_max"].isna()
    return frame


def robust_block(frame: pd.DataFrame, features: tuple[str, ...]) -> np.ndarray:
    values = frame.loc[:, features].to_numpy(float)
    values = RobustScaler(quantile_range=(10.0, 90.0)).fit_transform(values)
    values = np.clip(values, -8.0, 8.0)
    return values / math.sqrt(values.shape[1])


def knn_affinity(values: np.ndarray, k: int = K_GRAPH) -> sparse.csr_matrix:
    n = len(values)
    kk = min(k + 1, n)
    distances, indices = NearestNeighbors(
        n_neighbors=kk, metric="euclidean", n_jobs=-1
    ).fit(values).kneighbors(values)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    local_scale = np.maximum(distances[:, -1], 1e-10)
    rows = np.repeat(np.arange(n), indices.shape[1])
    cols = indices.ravel()
    denominator = np.sqrt(local_scale[rows] * local_scale[cols])
    weights = np.exp(-np.square(distances.ravel() / np.maximum(denominator, 1e-10)))
    graph = sparse.coo_matrix((weights, (rows, cols)), shape=(n, n)).tocsr()
    graph = graph.maximum(graph.T)
    graph.setdiag(0.0)
    graph.eliminate_zeros()
    return graph


def diffusion_coordinates(graph: sparse.csr_matrix, n_components: int = N_DIFFUSION) -> np.ndarray:
    degree = np.asarray(graph.sum(axis=1)).ravel()
    invsqrt = sparse.diags(1.0 / np.sqrt(np.maximum(degree, 1e-12)))
    normalized = invsqrt @ graph @ invsqrt
    k = min(n_components + 1, graph.shape[0] - 1)
    eigenvalues, eigenvectors = sparse.linalg.eigsh(
        normalized, k=k, which="LA", tol=1e-6
    )
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    coordinates = eigenvectors[:, 1:k] * eigenvalues[None, 1:k]
    for column in range(coordinates.shape[1]):
        pivot = int(np.argmax(np.abs(coordinates[:, column])))
        if coordinates[pivot, column] < 0:
            coordinates[:, column] *= -1
    return coordinates


def build_views(
    frame: pd.DataFrame,
) -> tuple[dict[str, np.ndarray], dict[str, sparse.csr_matrix], dict[str, np.ndarray]]:
    structure = robust_block(frame, STRUCTURE_FEATURES)
    electronic_base = robust_block(frame, ELECTRONIC_BASE_FEATURES)
    electronic_refined = robust_block(frame, ELECTRONIC_REFINED_FEATURES)
    values = {
        "structure S1": structure,
        "electronic E0": electronic_base,
        "electronic E2": electronic_refined,
        "joint S1+E0": np.hstack([structure / math.sqrt(2.0), electronic_base / math.sqrt(2.0)]),
        "joint S1+E2": np.hstack([structure / math.sqrt(2.0), electronic_refined / math.sqrt(2.0)]),
    }
    graphs = {name: knn_affinity(block) for name, block in values.items()}
    coordinates = {name: diffusion_coordinates(graph) for name, graph in graphs.items()}
    visual_names = ("structure S1", "electronic E2", "joint S1+E2")
    visual = {
        name: UMAP(
            n_neighbors=K_GRAPH,
            min_dist=0.18,
            metric="euclidean",
            random_state=SEED,
            n_jobs=1,
        ).fit_transform(coordinates[name])
        for name in visual_names
    }
    return coordinates, graphs, visual


def score_from_seed_distance(
    coordinates: np.ndarray,
    frame: pd.DataFrame,
    excluded_formula: str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seed = frame["seed_formula"].to_numpy(bool, copy=True)
    if excluded_formula is not None:
        seed &= frame["canon"].astype(str).to_numpy() != str(excluded_formula)
    seed_indices = np.flatnonzero(seed)
    if len(seed_indices) == 0:
        raise ValueError("No seed formulas remain")
    distance = cdist(coordinates, coordinates[seed_indices], metric="euclidean")
    nearest_column = np.argmin(distance, axis=1)
    nearest_distance = distance[np.arange(len(frame)), nearest_column]
    nearest_index = seed_indices[nearest_column]
    nearest_formula = frame.iloc[nearest_index]["canon"].astype(str).to_numpy()
    score = pd.Series(-nearest_distance).rank(pct=True, method="average").to_numpy()
    return score, nearest_distance, nearest_formula


def leave_one_formula_out(
    coordinates: dict[str, np.ndarray], frame: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict] = []
    formulas = sorted(frame.loc[frame["seed_formula"], "canon"].dropna().unique())
    for view, coords in coordinates.items():
        for held_formula in formulas:
            score, _, _ = score_from_seed_distance(coords, frame, held_formula)
            table = pd.DataFrame(
                {
                    "canon": frame["canon"].astype(str),
                    "score": score,
                    "training_seed": frame["seed_formula"]
                    & (frame["canon"].astype(str) != str(held_formula)),
                }
            )
            by_formula = table.groupby("canon", as_index=False).agg(
                score=("score", "max"), training_seed=("training_seed", "max")
            )
            eligible = by_formula.loc[~by_formula["training_seed"]].copy()
            eligible["retrieval_percentile"] = eligible["score"].rank(
                pct=True, method="average"
            )
            held = eligible.loc[eligible["canon"] == str(held_formula)]
            if len(held) != 1:
                continue
            percentile = float(held.iloc[0]["retrieval_percentile"])
            rows.append(
                {
                    "view": view,
                    "held_formula": held_formula,
                    "retrieval_percentile": percentile,
                    "recovered_top10pct": percentile >= 0.90,
                }
            )
    return pd.DataFrame(rows)


def add_screen_scores(frame: pd.DataFrame, joint_coordinates: np.ndarray) -> pd.DataFrame:
    out = frame.copy()
    analog, distance, nearest_formula = score_from_seed_distance(joint_coordinates, out)
    out["manifold_analog_percentile"] = analog
    out["distance_to_known_te_manifold"] = distance
    out["nearest_known_high_zt_formula"] = nearest_formula
    out["manifold_dual_screen_score"] = np.sqrt(
        out["manifold_analog_percentile"] * out["dual_score"].clip(0.0, 1.0)
    )
    candidate_pool = out["unknown_to_local_zt_table"] & ~out["seed_formula"]
    out["unknown_screen_percentile"] = np.nan
    out.loc[candidate_pool, "unknown_screen_percentile"] = out.loc[
        candidate_pool, "manifold_dual_screen_score"
    ].rank(pct=True, method="average")
    out["top_manifold_candidate"] = False
    selected = out.loc[candidate_pool, "manifold_dual_screen_score"].nlargest(
        min(N_CANDIDATES_PLOT, int(candidate_pool.sum()))
    ).index
    out.loc[selected, "top_manifold_candidate"] = True
    return out


def annotate_formulas(ax, subset: pd.DataFrame, x: str, y: str, n: int, fontsize: float = 8.0) -> None:
    labels = subset.drop_duplicates("canon").head(n)
    offsets = [(5, 5), (-5, 6), (5, -10), (-5, -11), (8, 2), (-8, 2)]
    for i, row in enumerate(labels.itertuples(index=False)):
        dx, dy = offsets[i % len(offsets)]
        ax.annotate(
            formula_label(row.canon),
            (getattr(row, x), getattr(row, y)),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=fontsize,
            ha="left" if dx > 0 else "right",
            va="bottom" if dy > 0 else "top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.25},
            zorder=9,
        )


def plot(frame: pd.DataFrame, validation: pd.DataFrame, summary: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 12.2), constrained_layout=True)
    map_specs = [
        (axes[0, 0], "structure_phi1", "structure_phi2", "Structure manifold: S1 geometry + composition"),
        (axes[0, 1], "electronic_phi1", "electronic_phi2", "Electronic manifold: E2 gap + dielectric + effective mass"),
        (axes[1, 0], "joint_phi1", "joint_phi2", "Balanced joint diffusion manifold: S1 × E2"),
    ]
    other = frame[~frame["seed_formula"] & ~frame["top_manifold_candidate"]]
    candidates = frame[frame["top_manifold_candidate"]].sort_values(
        "manifold_dual_screen_score", ascending=False
    )
    known = frame[frame["seed_formula"]].sort_values("external_zt_max", ascending=False)
    for ax, x, y, title in map_specs:
        ax.scatter(
            other[x], other[y], s=8, color=COLORS["other"], alpha=0.24,
            edgecolors="none", rasterized=True, zorder=1,
        )
        ax.scatter(
            candidates[x], candidates[y], s=48, marker="o", color=COLORS["candidate"],
            edgecolors="#222222", linewidths=0.45, alpha=0.95, zorder=5,
        )
        ax.scatter(
            known[x], known[y], s=88, marker="*", color=COLORS["known"],
            edgecolors="#111111", linewidths=0.65, zorder=7,
        )
        ax.set_xlabel("2D manifold layout 1")
        ax.set_ylabel("2D manifold layout 2")
        ax.set_title(title)
        ax.tick_params(labelsize=8)
    annotate_formulas(axes[1, 0], known, "joint_phi1", "joint_phi2", 8, 7.8)
    annotate_formulas(axes[1, 0], candidates, "joint_phi1", "joint_phi2", 6, 7.6)

    ax = axes[1, 1]
    order = [
        "structure S1",
        "electronic E0",
        "electronic E2",
        "joint S1+E0",
        "joint S1+E2",
    ]
    labels = ["Structure S1", "Electronic E0", "Electronic E2", "Joint S1+E0", "Joint S1+E2"]
    grouped = validation.groupby("view").agg(
        median_percentile=("retrieval_percentile", "median"),
        recall_top10=("recovered_top10pct", "mean"),
        n_formulas=("held_formula", "nunique"),
    ).loc[order]
    colors = [
        COLORS["structure"], COLORS["electronic"], COLORS["electronic"],
        COLORS["joint"], COLORS["joint"],
    ]
    bars = ax.barh(np.arange(len(order))[::-1], grouped["median_percentile"], color=colors, alpha=0.88)
    ax.axvline(0.5, color="#666666", lw=1.0, ls="--")
    ax.set_yticks(np.arange(len(order))[::-1], labels)
    ax.set_xlim(0.0, 1.04)
    ax.set_xlabel("held-out high-zT formula retrieval percentile")
    ax.set_title("Leave-one-formula-out manifold test\nfull-dimensional diffusion coordinates")
    for bar, (_, row) in zip(bars, grouped.iterrows()):
        ax.text(
            min(1.015, float(row.median_percentile) + 0.018),
            bar.get_y() + bar.get_height() / 2,
            f"median={row.median_percentile:.2f}; top-10% {int(round(row.recall_top10 * row.n_formulas))}/{int(row.n_formulas)}",
            va="center", fontsize=8.4,
        )
    ax.text(0.505, -0.65, "random median", fontsize=8.2, color="#555555", ha="left")

    legend = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["other"], markeredgecolor="none", markersize=6, label="other complete-case materials"),
        Line2D([0], [0], marker="*", linestyle="none", markerfacecolor=COLORS["known"], markeredgecolor="#111111", markersize=11, label=f"reported max zT≥1 formula match ({summary['n_seed_formulas']} formulas)"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["candidate"], markeredgecolor="#222222", markersize=8, label=f"top {summary['n_candidates_highlighted']} unknown analogues"),
    ]
    axes[0, 0].legend(handles=legend, loc="best", fontsize=8.2, framealpha=0.88)
    fig.suptitle(
        "Thermoelectric analogue screening on target-free structure–electronic manifolds\n"
        "Graph coordinates exclude PF, kL and zT; purple ranking combines joint-manifold analogy with existing dual-channel consistency",
        fontsize=14,
    )
    fig.savefig(FIGURE_OUT, dpi=230, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    module = load_descriptor_module()
    frame = prepare_data(module)
    coordinates, graphs, visual = build_views(frame)
    validation = leave_one_formula_out(coordinates, frame)
    frame = add_screen_scores(frame, coordinates["joint S1+E2"])

    coordinate_map = {
        "structure S1": "structure",
        "electronic E2": "electronic",
        "joint S1+E2": "joint",
    }
    for view, prefix in coordinate_map.items():
        frame[f"{prefix}_phi1"] = visual[view][:, 0]
        frame[f"{prefix}_phi2"] = visual[view][:, 1]

    candidates = frame[
        frame["unknown_to_local_zt_table"] & ~frame["seed_formula"]
    ].sort_values("manifold_dual_screen_score", ascending=False)
    candidate_columns = [
        "row_id", "formula", "canon", "chemical_system", "preferred_carrier",
        "nearest_known_high_zt_formula", "manifold_analog_percentile",
        "structure_score_percentile", "electronic_score_percentile", "dual_score",
        "manifold_dual_screen_score", "unknown_screen_percentile",
    ]
    candidates[candidate_columns].head(100).to_csv(CANDIDATES_OUT, index=False)

    validation.to_csv(VALIDATION_OUT, index=False)
    point_columns = [
        "row_id", "formula", "canon", "chemical_system", "preferred_carrier",
        "seed_formula", "external_zt_max", "unknown_to_local_zt_table",
        "nearest_known_high_zt_formula", "manifold_analog_percentile",
        "distance_to_known_te_manifold", "structure_score_percentile",
        "electronic_score_percentile", "dual_score", "manifold_dual_screen_score",
        "unknown_screen_percentile", "top_manifold_candidate",
        "structure_phi1", "structure_phi2", "electronic_phi1", "electronic_phi2",
        "joint_phi1", "joint_phi2",
    ]
    frame[point_columns].to_csv(POINTS_OUT, index=False)

    metric = validation.groupby("view").agg(
        median_retrieval_percentile=("retrieval_percentile", "median"),
        mean_retrieval_percentile=("retrieval_percentile", "mean"),
        top10_recall=("recovered_top10pct", "mean"),
        n_seed_formulas=("held_formula", "nunique"),
    )
    components = {name: int(graph.shape[0]) for name, graph in graphs.items()}
    summary = {
        "no_new_first_principles_or_transport_calculation": True,
        "manifold_construction_is_target_free": True,
        "n_complete_case_materials": int(len(frame)),
        "n_seed_jids": int(frame["seed_formula"].sum()),
        "n_seed_formulas": int(frame.loc[frame["seed_formula"], "canon"].nunique()),
        "n_unknown_formula_table_materials": int(frame["unknown_to_local_zt_table"].sum()),
        "n_candidates_highlighted": int(frame["top_manifold_candidate"].sum()),
        "graph_k": K_GRAPH,
        "diffusion_dimensions_for_scoring": N_DIFFUSION,
        "structure_features": list(STRUCTURE_FEATURES),
        "electronic_features": list(ELECTRONIC_REFINED_FEATURES),
        "view_metrics": metric.to_dict(orient="index"),
        "graph_node_counts": components,
        "screening_definition": (
            "geometric mean of full-dimensional joint-manifold seed-proximity percentile "
            "and the pre-existing cross-validated dual-channel score"
        ),
        "warning": (
            "high-zT seeds are reduced-formula matches; polymorph, doping, dimensionality, "
            "microstructure and measurement conditions may differ"
        ),
    }
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    plot(frame, validation, summary)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nTop candidates")
    print(candidates[candidate_columns].head(30).to_string(index=False))


if __name__ == "__main__":
    main()
