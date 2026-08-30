"""Global, reference-free structure/electronic common-space visualization.

Every coordinate is constructed from all materials before PF, experimental
kappa_L, or external zT labels are inspected.  The cyan/purple/orange marks are
post-hoc overlays only.  No material is used as an origin or reference point.
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
from scipy.linalg import orthogonal_procrustes


HERE = Path(__file__).resolve().parent
BASE_SCRIPT = HERE / "run_rich_descriptor_and_space.py"
OUTPUT_DIR = HERE / "outputs"
FIGURE_DIR = HERE / "figures"
POINTS_OUT = OUTPUT_DIR / "global_common_space_points.csv"
SUMMARY_OUT = OUTPUT_DIR / "global_common_space_summary.json"
FIGURE_OUT = FIGURE_DIR / "global_structure_electronic_common_space.png"
PDF_OUT = FIGURE_DIR / "global_structure_electronic_common_space.pdf"

GREY = "#c9cdd2"
CYAN = "#00cbd8"
CYAN_EDGE = "#008c96"
PURPLE = "#8e44ad"
DARK_PURPLE = "#54226d"
ORANGE = "#e07a16"


def load_base_module():
    spec = importlib.util.spec_from_file_location("rich_descriptor_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(BASE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def centre_and_scale(coordinates: np.ndarray) -> np.ndarray:
    result = coordinates - coordinates.mean(axis=0, keepdims=True)
    scale = np.sqrt(np.mean(np.sum(result**2, axis=1)))
    return result / max(scale, 1e-12)


def align_to(reference: np.ndarray, moving: np.ndarray) -> np.ndarray:
    reference = centre_and_scale(reference)
    moving = centre_and_scale(moving)
    rotation, _ = orthogonal_procrustes(moving, reference)
    return moving @ rotation


def add_global_coordinates(base, frame: pd.DataFrame, distances: dict[str, np.ndarray]):
    coordinate_sets: dict[str, np.ndarray] = {}
    preservation: dict[str, float] = {}
    for name in ["structure", "pure_electronic", "pure_joint_and"]:
        coordinate_sets[name], preservation[name] = base.pcoa_mds(distances[name])
    coordinate_sets["structure"] = centre_and_scale(coordinate_sets["structure"])
    coordinate_sets["pure_electronic"] = align_to(
        coordinate_sets["structure"], coordinate_sets["pure_electronic"]
    )
    coordinate_sets["pure_joint_and"] = align_to(
        coordinate_sets["structure"], coordinate_sets["pure_joint_and"]
    )
    result = frame.copy()
    for name, coordinates in coordinate_sets.items():
        result[f"{name}_x"] = coordinates[:, 0]
        result[f"{name}_y"] = coordinates[:, 1]
    return result, preservation


def draw_panel(ax, frame: pd.DataFrame, prefix: str, title: str, preservation: float) -> None:
    x = f"{prefix}_x"
    y = f"{prefix}_y"
    grey = frame[~frame["known_high_zt"] & ~frame["dual_screening_hit"]]
    high = frame[frame["known_high_zt"]]
    screen = frame[frame["dual_screening_hit"]]
    reported_low = frame[frame["reported_low_zt_hit"]]
    overlap = frame[frame["known_high_zt"] & frame["dual_good"]]

    ax.scatter(grey[x], grey[y], s=24, c=GREY, alpha=0.58, edgecolors="none", zorder=1)
    ax.scatter(screen[x], screen[y], s=82, c=PURPLE, edgecolors="white", linewidths=0.8, zorder=5)
    if len(reported_low):
        ax.scatter(reported_low[x], reported_low[y], s=43, marker="x", c=ORANGE, linewidths=1.25, zorder=7)
    ax.scatter(high[x], high[y], s=145, marker="*", c=CYAN, edgecolors=CYAN_EDGE, linewidths=0.8, zorder=6)
    if len(overlap):
        ax.scatter(overlap[x], overlap[y], s=205, facecolors="none", edgecolors=DARK_PURPLE, linewidths=1.7, zorder=8)

    for _, row in screen.iterrows():
        dy = -12 if row[y] > frame[y].quantile(0.9) else 5
        ax.annotate(row["formula"], (row[x], row[y]), xytext=(4, dy), textcoords="offset points", fontsize=7.8, color=DARK_PURPLE)
    for _, row in high.sort_values("external_zt_max", ascending=False).head(4).iterrows():
        ax.annotate(row["formula"], (row[x], row[y]), xytext=(4, -11), textcoords="offset points", fontsize=7.3, color=CYAN_EDGE)

    ax.set_xlabel("global map coordinate 1 (dimensionless)")
    ax.set_ylabel("global map coordinate 2 (dimensionless)")
    ax.set_title(f"{title}\npairwise-rank preservation = {preservation:.2f}")
    ax.grid(color="#eeeeee", lw=0.6, zorder=0)


def plot(frames: dict[str, pd.DataFrame], preservations: dict[str, dict[str, float]]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17.2, 10.8))
    panels = [
        ("structure", "global structure space"),
        ("pure_electronic", "global electronic-structure space"),
        ("pure_joint_and", "global strict-AND common space"),
    ]
    for row, carrier in enumerate(["n", "p"]):
        for column, (prefix, title) in enumerate(panels):
            draw_panel(
                axes[row, column],
                frames[carrier],
                prefix,
                f"{carrier}-type | {title}",
                preservations[carrier][prefix],
            )

    handles = [
        Line2D([], [], marker="o", ls="", color=GREY, label="other material"),
        Line2D([], [], marker="*", ls="", markersize=13, markerfacecolor=CYAN, markeredgecolor=CYAN_EDGE, label="formula-matched experimental zT >= 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor=PURPLE, markeredgecolor="white", label="PF-low-kL dual screening hit"),
        Line2D([], [], marker="x", ls="", color=ORANGE, label="external formula-level zT report < 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=DARK_PURPLE, markeredgewidth=1.6, label="known high-zT also passes dual screen"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.008))
    fig.suptitle(
        "Reference-free global representation of all materials\n"
        "coordinates are built before zT, PF, and experimental kappa_L are overlaid",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.085,
        "Strict-AND distance = max(global structure neighbour-rank, global electronic-structure neighbour-rank); no seed and no 1/2 block weight",
        ha="center",
        fontsize=9.5,
    )
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.17, top=0.86, wspace=0.22, hspace=0.34)
    fig.savefig(FIGURE_OUT, dpi=240, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    base = load_base_module()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    views = pd.read_parquet(base.VIEWS_PATH).reset_index(drop=True)
    soap = np.load(base.SOAP_PATH)
    raw = json.load(open(base.RAW_PATH))
    raw_by_jid = {str(record["jid"]): record for record in raw}
    dual = pd.read_csv(base.DUAL_PATH)
    zt_by_formula = base.external_zt_by_formula()

    frames: dict[str, pd.DataFrame] = {}
    preservations: dict[str, dict[str, float]] = {}
    results = {}
    for carrier in ["n", "p"]:
        frame, _, distances = base.build_carrier_frame(
            carrier, views, soap, raw_by_jid, zt_by_formula, dual
        )
        frame, preservation = add_global_coordinates(base, frame, distances)
        frames[carrier] = frame
        preservations[carrier] = preservation
        labels = frame["known_high_zt"].to_numpy(bool)
        formulas = frame["canon"].astype(str).to_numpy()
        results[carrier] = {
            "n_materials": int(len(frame)),
            "n_high_zt_formula_marks": int(labels.sum()),
            "n_dual_screen_hits_not_high_zt": int(frame["dual_screening_hit"].sum()),
            "n_known_high_zt_also_dual_screen": int((frame["known_high_zt"] & frame["dual_good"]).sum()),
            "lofo_auc_structure": base.lofo_auc(labels, formulas, labels, distances["structure"]),
            "lofo_auc_electronic_structure": base.lofo_auc(labels, formulas, labels, distances["pure_electronic"]),
            "lofo_auc_strict_and_common": base.lofo_auc(labels, formulas, labels, distances["pure_joint_and"]),
            "map_pairwise_rank_preservation": preservation,
        }

    points = pd.concat(frames.values(), ignore_index=True)
    point_columns = [
        "carrier", "jid", "formula", "canon", "known_high_zt", "external_zt_max",
        "dual_good", "dual_screening_hit", "reported_low_zt_hit", "unlabelled_dual_hit",
        "PF_percentile", "low_kL_percentile",
        "structure_x", "structure_y", "pure_electronic_x", "pure_electronic_y",
        "pure_joint_and_x", "pure_joint_and_y",
    ]
    points[point_columns].to_csv(POINTS_OUT, index=False)
    summary = {
        "reference_free": True,
        "labels_used_in_coordinates": [],
        "joint_distance": "max(symmetrized structure neighbour-rank, symmetrized pure-electronic neighbour-rank)",
        "electronic_space_excludes": ["S", "sigma", "kappa_e", "PF", "zT"],
        "results": results,
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    plot(frames, preservations)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Saved {FIGURE_OUT}")


if __name__ == "__main__":
    main()
