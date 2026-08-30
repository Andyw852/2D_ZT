"""Build an interpretable structure/electronic strict-AND space from local data.

No DFT, BTE, phonon, or new transport calculation is performed.  The script
only reuses the existing JARVIS/StarryData2 tables and the already-computed
geometry-only SOAP matrix.

Two electronic views are deliberately kept separate:

* band view: gap, effective masses, dielectric response, and spillage;
* transport fingerprint: fixed-condition S, sigma, and kappa_e.

The transport fingerprint is useful for analogue visualization, but it is not
an independent validation of PF because PF is algebraically related to S and
sigma.  PF, experimental kappa_L, and external zT are never used to construct
the coordinates.  They are only used after construction for colours/labels.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pymatgen.core import Composition, Element
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances, roc_auc_score
from sklearn.metrics.pairwise import cosine_distances
from sklearn.preprocessing import RobustScaler


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[1]
ATLAS = WORKSPACE / "jarvis_2d_te_atlas"

VIEWS_PATH = ATLAS / "features/kl_verify/kl_views.parquet"
SOAP_PATH = ATLAS / "data/processed/kl_soap_geo.npy"
RAW_PATH = ATLAS / "data/raw/external/jarvis_kl/jdft_3d-8-18-2021.json"
DUAL_PATH = ATLAS / "data/processed/pf_kL_dual_channel_intersection.csv"
ZT_PATH = ATLAS / "data/raw/external/starrydata2/starrydata_zt_sane.csv"

OUTPUT_DIR = HERE / "outputs"
FIGURE_DIR = HERE / "figures"
POINTS_OUT = OUTPUT_DIR / "rich_descriptor_points.csv"
CANDIDATES_OUT = OUTPUT_DIR / "rich_descriptor_candidates.csv"
SUMMARY_OUT = OUTPUT_DIR / "rich_descriptor_summary.json"
SIMILARITY_FIGURE = FIGURE_DIR / "rich_descriptor_similarity_planes.png"
SIMILARITY_PDF = FIGURE_DIR / "rich_descriptor_similarity_planes.pdf"
AND_FIGURE = FIGURE_DIR / "rich_descriptor_strict_and_map.png"
AND_PDF = FIGURE_DIR / "rich_descriptor_strict_and_map.pdf"

RANDOM_STATE = 20260830
HIGH_ZT_THRESHOLD = 1.0
ANALOGUE_THRESHOLD = 0.75

GREY = "#c9cdd2"
CYAN = "#00cbd8"
CYAN_EDGE = "#008c96"
PURPLE = "#8e44ad"
DARK_PURPLE = "#54226d"
ORANGE = "#e07a16"


def canon(formula: object) -> str | None:
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return None


def finite(value: object, *, positive: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return np.nan
    if not np.isfinite(result) or result <= -99998:
        return np.nan
    if positive and result <= 0:
        return np.nan
    return result


def weighted_stats(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    good = np.isfinite(values) & np.isfinite(weights)
    if not good.any():
        return np.nan, np.nan, np.nan
    values = values[good]
    weights = weights[good]
    weights = weights / weights.sum()
    mean = float(np.sum(values * weights))
    std = float(np.sqrt(np.sum(weights * (values - mean) ** 2)))
    return mean, std, float(values.max() - values.min())


def composition_features(formula: str) -> dict[str, float]:
    comp = Composition(formula).fractional_composition
    items = list(comp.items())
    weights = np.asarray([float(amount) for _, amount in items], dtype=float)
    weights /= weights.sum()
    elements = [Element(str(element)) for element, _ in items]

    property_values: dict[str, np.ndarray] = {
        "atomic_mass": np.asarray([finite(element.atomic_mass) for element in elements]),
        "atomic_number": np.asarray([float(element.Z) for element in elements]),
        "electronegativity": np.asarray([finite(element.X) for element in elements]),
        "atomic_radius": np.asarray([finite(element.atomic_radius) for element in elements]),
        "mendeleev": np.asarray([finite(element.mendeleev_no) for element in elements]),
    }
    result: dict[str, float] = {
        "n_elements": float(len(elements)),
        "composition_entropy": float(-np.sum(weights * np.log(weights + 1e-15))),
        "largest_element_fraction": float(weights.max()),
    }
    for name, values in property_values.items():
        mean, std, value_range = weighted_stats(values, weights)
        result[f"{name}_mean"] = mean
        result[f"{name}_std"] = std
        result[f"{name}_range"] = value_range
    return result


def lattice_features(record: dict) -> dict[str, float]:
    atoms = record.get("atoms") or {}
    lattice = np.asarray(atoms.get("lattice_mat", np.full((3, 3), np.nan)), dtype=float)
    nat = len(atoms.get("elements", []))
    lengths = np.linalg.norm(lattice, axis=1) if lattice.shape == (3, 3) else np.full(3, np.nan)
    volume = abs(float(np.linalg.det(lattice))) if lattice.shape == (3, 3) else np.nan
    volume_per_atom = volume / nat if nat and np.isfinite(volume) and volume > 0 else np.nan

    def angle(u: np.ndarray, v: np.ndarray) -> float:
        denom = np.linalg.norm(u) * np.linalg.norm(v)
        if denom <= 0 or not np.isfinite(denom):
            return np.nan
        return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / denom, -1.0, 1.0))))

    angles = (
        np.asarray([angle(lattice[1], lattice[2]), angle(lattice[0], lattice[2]), angle(lattice[0], lattice[1])])
        if lattice.shape == (3, 3)
        else np.full(3, np.nan)
    )
    b_modulus = finite(record.get("bulk_modulus_kv"), positive=True)
    g_modulus = finite(record.get("shear_modulus_gv"), positive=True)
    density = finite(record.get("density"), positive=True)
    feature = {
        "log_density": math.log10(density) if np.isfinite(density) else np.nan,
        "log_volume_per_atom": math.log10(volume_per_atom) if np.isfinite(volume_per_atom) else np.nan,
        "log_n_atoms": math.log10(nat) if nat > 0 else np.nan,
        "length_cv": float(np.nanstd(lengths) / np.nanmean(lengths)) if np.isfinite(lengths).any() else np.nan,
        "angle_deviation": float(np.nanmean(np.abs(angles - 90.0))) if np.isfinite(angles).any() else np.nan,
        "log_bulk_modulus": math.log10(b_modulus) if np.isfinite(b_modulus) else np.nan,
        "log_shear_modulus": math.log10(g_modulus) if np.isfinite(g_modulus) else np.nan,
        "log_B_over_G": math.log10(b_modulus / g_modulus) if np.isfinite(b_modulus) and np.isfinite(g_modulus) else np.nan,
        "poisson": finite(record.get("poisson")),
    }
    spg = finite(record.get("spg_number"), positive=True)
    crystal_system = "unknown"
    if np.isfinite(spg):
        if spg <= 2:
            crystal_system = "triclinic"
        elif spg <= 15:
            crystal_system = "monoclinic"
        elif spg <= 74:
            crystal_system = "orthorhombic"
        elif spg <= 142:
            crystal_system = "tetragonal"
        elif spg <= 167:
            crystal_system = "trigonal"
        elif spg <= 194:
            crystal_system = "hexagonal"
        else:
            crystal_system = "cubic"
    for system in ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"]:
        feature[f"crystal_{system}"] = float(crystal_system == system)
    return feature


def band_features(record: dict, carrier: str) -> dict[str, float]:
    eg_opt = finite(record.get("optb88vdw_bandgap"))
    eg_mbj = finite(record.get("mbj_bandgap"))
    m_elec = finite(record.get("avg_elec_mass"), positive=True)
    m_hole = finite(record.get("avg_hole_mass"), positive=True)
    carrier_mass = m_elec if carrier == "n" else m_hole
    opposite_mass = m_hole if carrier == "n" else m_elec
    eps = np.asarray([finite(record.get(key), positive=True) for key in ["epsx", "epsy", "epsz"]])
    good_eps = eps[np.isfinite(eps)]
    eps_geo = float(np.exp(np.mean(np.log(good_eps)))) if len(good_eps) else np.nan
    eps_anisotropy = float(good_eps.max() / good_eps.min()) if len(good_eps) >= 2 else np.nan
    return {
        "Eg_opt": eg_opt,
        "Eg_mbj": eg_mbj,
        "gap_correction_mbj_minus_opt": eg_mbj - eg_opt if np.isfinite(eg_mbj) and np.isfinite(eg_opt) else np.nan,
        "log_carrier_mass": math.log10(carrier_mass) if np.isfinite(carrier_mass) else np.nan,
        "log_opposite_mass": math.log10(opposite_mass) if np.isfinite(opposite_mass) else np.nan,
        "log_mass_ratio": math.log10(carrier_mass / opposite_mass) if np.isfinite(carrier_mass) and np.isfinite(opposite_mass) else np.nan,
        "log_dielectric_geo": math.log10(eps_geo) if np.isfinite(eps_geo) else np.nan,
        "log_dielectric_anisotropy": math.log10(eps_anisotropy) if np.isfinite(eps_anisotropy) else np.nan,
        "log1p_spillage": math.log10(1.0 + finite(record.get("spillage"), positive=True)) if np.isfinite(finite(record.get("spillage"), positive=True)) else np.nan,
    }


def transport_features(record: dict, carrier: str) -> dict[str, float]:
    seebeck = finite(record.get(f"{carrier}-Seebeck"))
    sigma = finite(record.get(f"{carrier}cond"), positive=True)
    kappa_e = finite(record.get(f"{carrier}kappa"), positive=True)
    expected_sign = -1.0 if carrier == "n" else 1.0
    return {
        "abs_Seebeck": abs(seebeck) if np.isfinite(seebeck) else np.nan,
        "Seebeck_sign_expected": float(np.sign(seebeck) == expected_sign) if np.isfinite(seebeck) else np.nan,
        "log_sigma": math.log10(sigma) if np.isfinite(sigma) else np.nan,
        "log_kappa_e": math.log10(kappa_e) if np.isfinite(kappa_e) else np.nan,
        "log_kappa_e_over_sigma": math.log10(kappa_e / sigma) if np.isfinite(kappa_e) and np.isfinite(sigma) else np.nan,
    }


def numeric_distance(frame: pd.DataFrame) -> np.ndarray:
    values = frame.to_numpy(dtype=float)
    transformed = SimpleImputer(strategy="median", add_indicator=True).fit_transform(values)
    transformed = RobustScaler(quantile_range=(10.0, 90.0)).fit_transform(transformed)
    return pairwise_distances(transformed, metric="euclidean") / math.sqrt(max(1, transformed.shape[1]))


def row_rank_matrix(distance: np.ndarray) -> np.ndarray:
    n = len(distance)
    order = np.argsort(distance, axis=1, kind="mergesort")
    ranks = np.empty((n, n), dtype=float)
    ranks[np.arange(n)[:, None], order] = np.arange(n, dtype=float)[None, :]
    ranks /= max(1, n - 1)
    np.fill_diagonal(ranks, 0.0)
    return ranks


def pcoa_mds(distance: np.ndarray) -> tuple[np.ndarray, float]:
    model = MDS(
        n_components=2,
        metric=True,
        dissimilarity="precomputed",
        random_state=RANDOM_STATE,
        n_init=16,
        max_iter=1000,
        eps=1e-8,
        normalized_stress="auto",
        init="classical_mds",
    )
    coordinates = model.fit_transform(distance)
    embedded = pairwise_distances(coordinates)
    triangle = np.triu_indices(len(distance), 1)
    rho = float(spearmanr(distance[triangle], embedded[triangle]).statistic)
    return coordinates, rho


def external_zt_by_formula() -> pd.DataFrame:
    zt = pd.read_csv(ZT_PATH)
    zt = zt[(zt["T_max"] >= 200.0) & (zt["T_max"] <= 1500.0)].copy()
    zt["canon"] = zt["composition"].map(canon)
    zt = zt.dropna(subset=["canon", "zt_max"])
    return zt.groupby("canon", as_index=False).agg(
        external_zt_max=("zt_max", "max"),
        external_zt_samples=("sample_id", "nunique"),
    )


def nearest_seed_scores(
    formula: np.ndarray,
    seed: np.ndarray,
    structure_rank: np.ndarray,
    electronic_rank: np.ndarray,
) -> pd.DataFrame:
    seed_index = np.flatnonzero(seed)
    rows = []
    for index in range(len(formula)):
        eligible = seed_index
        if seed[index]:
            eligible = seed_index[formula[seed_index] != formula[index]]
        if len(eligible) == 0:
            rows.append((np.nan, np.nan, np.nan, -1))
            continue
        pair_worst = np.maximum(structure_rank[index, eligible], electronic_rank[index, eligible])
        selected = int(eligible[int(np.argmin(pair_worst))])
        s_sim = 1.0 - float(structure_rank[index, selected])
        e_sim = 1.0 - float(electronic_rank[index, selected])
        rows.append((s_sim, e_sim, min(s_sim, e_sim), selected))
    return pd.DataFrame(rows, columns=["structure_similarity", "electronic_similarity", "and_similarity", "selected_seed_index"])


def lofo_auc(labels: np.ndarray, formula: np.ndarray, seed: np.ndarray, distance_rank: np.ndarray) -> float:
    scores = []
    for index in range(len(labels)):
        eligible = np.flatnonzero(seed)
        if labels[index]:
            eligible = eligible[formula[eligible] != formula[index]]
        scores.append(1.0 - float(distance_rank[index, eligible].min()) if len(eligible) else np.nan)
    good = np.isfinite(scores)
    return float(roc_auc_score(labels[good], np.asarray(scores)[good]))


def build_carrier_frame(
    carrier: str,
    views: pd.DataFrame,
    soap_all: np.ndarray,
    raw_by_jid: dict[str, dict],
    zt_by_formula: pd.DataFrame,
    dual: pd.DataFrame,
) -> tuple[pd.DataFrame, dict, dict[str, np.ndarray]]:
    frame = dual[dual["carrier"] == carrier].copy().reset_index(drop=True)
    frame["canon"] = frame["formula"].map(canon)
    frame = frame.merge(zt_by_formula, on="canon", how="left", validate="many_to_one")
    frame["known_high_zt"] = frame["external_zt_max"].fillna(-np.inf) >= HIGH_ZT_THRESHOLD
    frame["dual_good"] = frame["top20_intersection"].astype(bool) | frame["pareto"].astype(bool)
    # Some non-high-zT screening hits already have an external formula-level
    # report below one; others are truly unlabelled.  Keep the distinction.
    frame["dual_screening_hit"] = frame["dual_good"] & ~frame["known_high_zt"]
    frame["reported_low_zt_hit"] = frame["dual_screening_hit"] & frame["external_zt_max"].notna()
    frame["unlabelled_dual_hit"] = frame["dual_screening_hit"] & frame["external_zt_max"].isna()
    frame["dual_desirability"] = np.sqrt(frame["PF_percentile"] * frame["low_kL_percentile"])

    view_index = {jid: index for index, jid in enumerate(views["jid"].astype(str))}
    source_index = np.asarray([view_index[jid] for jid in frame["jid"].astype(str)], dtype=int)
    soap = np.asarray(soap_all[source_index], dtype=float)
    soap_norm = np.linalg.norm(soap, axis=1, keepdims=True)
    soap = soap / np.where(soap_norm > 0, soap_norm, 1.0)

    composition_rows = [composition_features(formula) for formula in frame["formula"]]
    lattice_rows = [lattice_features(raw_by_jid[jid]) for jid in frame["jid"]]
    band_rows = [band_features(raw_by_jid[jid], carrier) for jid in frame["jid"]]
    transport_rows = [transport_features(raw_by_jid[jid], carrier) for jid in frame["jid"]]

    composition_lattice = pd.DataFrame(
        [{**composition, **lattice} for composition, lattice in zip(composition_rows, lattice_rows)]
    )
    band = pd.DataFrame(band_rows)
    transport = pd.DataFrame(transport_rows)

    geometry_rank = row_rank_matrix(cosine_distances(soap))
    lattice_rank = row_rank_matrix(numeric_distance(composition_lattice))
    band_rank = row_rank_matrix(numeric_distance(band))
    transport_rank = row_rank_matrix(numeric_distance(transport))

    # Strict conjunction: a poor rank in either sub-view cannot be compensated
    # by a good rank in the other.  No 1/2 or fitted block weight is used.
    structure_rank = np.maximum(geometry_rank, lattice_rank)
    rich_electronic_rank = np.maximum(band_rank, transport_rank)
    pure = nearest_seed_scores(
        frame["canon"].astype(str).to_numpy(),
        frame["known_high_zt"].to_numpy(bool),
        structure_rank,
        band_rank,
    ).add_prefix("pure_")
    rich = nearest_seed_scores(
        frame["canon"].astype(str).to_numpy(),
        frame["known_high_zt"].to_numpy(bool),
        structure_rank,
        rich_electronic_rank,
    ).add_prefix("rich_")
    frame = pd.concat([frame, pure, rich], axis=1)

    seed_formula = frame["canon"].astype(str).to_numpy()
    selected = frame["rich_selected_seed_index"].fillna(-1).astype(int).to_numpy()
    frame["selected_high_zt_seed"] = [seed_formula[index] if index >= 0 else "" for index in selected]
    frame["strict_analogue"] = frame["dual_screening_hit"] & (
        frame["rich_structure_similarity"] >= ANALOGUE_THRESHOLD
    ) & (frame["rich_electronic_similarity"] >= ANALOGUE_THRESHOLD)
    frame["geometry_similarity_to_seed"] = [
        1.0 - float(geometry_rank[row, index]) if index >= 0 else np.nan
        for row, index in enumerate(selected)
    ]
    frame["composition_lattice_similarity_to_seed"] = [
        1.0 - float(lattice_rank[row, index]) if index >= 0 else np.nan
        for row, index in enumerate(selected)
    ]
    frame["band_similarity_to_seed"] = [
        1.0 - float(band_rank[row, index]) if index >= 0 else np.nan
        for row, index in enumerate(selected)
    ]
    frame["transport_similarity_to_seed"] = [
        1.0 - float(transport_rank[row, index]) if index >= 0 else np.nan
        for row, index in enumerate(selected)
    ]
    frame["external_validation_status"] = np.select(
        [frame["known_high_zt"], frame["reported_low_zt_hit"], frame["unlabelled_dual_hit"]],
        ["known high-zT", "reported formula-level zT < 1", "unlabelled"],
        default="not a dual hit",
    )

    global_distances = {
        "structure": np.maximum(structure_rank, structure_rank.T),
        "pure_electronic": np.maximum(band_rank, band_rank.T),
        "rich_electronic": np.maximum(rich_electronic_rank, rich_electronic_rank.T),
    }
    global_distances["pure_joint_and"] = np.maximum(
        global_distances["structure"], global_distances["pure_electronic"]
    )
    global_distances["rich_joint_and"] = np.maximum(
        global_distances["structure"], global_distances["rich_electronic"]
    )
    for distance in global_distances.values():
        np.fill_diagonal(distance, 0.0)
    coordinates, preservation = pcoa_mds(global_distances["rich_joint_and"])
    frame["and_map_1"] = coordinates[:, 0]
    frame["and_map_2"] = coordinates[:, 1]

    labels = frame["known_high_zt"].to_numpy(bool)
    formulas = frame["canon"].astype(str).to_numpy()
    summary = {
        "carrier": carrier,
        "n_materials": int(len(frame)),
        "n_high_zt_formula_seeds": int(labels.sum()),
        "n_non_high_zt_dual_hits": int(frame["dual_screening_hit"].sum()),
        "n_reported_low_zt_hits": int(frame["reported_low_zt_hit"].sum()),
        "n_unlabelled_dual_hits": int(frame["unlabelled_dual_hit"].sum()),
        "n_strict_analogues": int(frame["strict_analogue"].sum()),
        "lofo_auc_structure_strict": lofo_auc(labels, formulas, labels, structure_rank),
        "lofo_auc_pure_band": lofo_auc(labels, formulas, labels, band_rank),
        "lofo_auc_transport_fingerprint": lofo_auc(labels, formulas, labels, transport_rank),
        "lofo_auc_rich_electronic_strict": lofo_auc(labels, formulas, labels, rich_electronic_rank),
        "lofo_auc_full_strict_and": lofo_auc(labels, formulas, labels, global_distances["rich_joint_and"]),
        "strict_and_map_pairwise_rank_preservation": preservation,
        "descriptor_counts": {
            "geometry_soap": int(soap.shape[1]),
            "composition_lattice": int(composition_lattice.shape[1]),
            "band": int(band.shape[1]),
            "transport_fingerprint": int(transport.shape[1]),
        },
    }
    return frame, summary, global_distances


def scatter_similarity(ax, frame: pd.DataFrame, prefix: str, carrier: str, title: str) -> None:
    grey = frame[~frame["known_high_zt"] & ~frame["dual_screening_hit"]]
    seed = frame[frame["known_high_zt"]]
    candidate = frame[frame["dual_screening_hit"]]
    reported_low = frame[frame["reported_low_zt_hit"]]
    strict = frame[frame["strict_analogue"]]
    x = f"{prefix}_structure_similarity"
    y = f"{prefix}_electronic_similarity"
    ax.scatter(grey[x], grey[y], s=18, c=GREY, alpha=0.55, edgecolors="none", zorder=1)
    ax.scatter(seed[x], seed[y], s=125, marker="*", facecolors=CYAN, edgecolors=CYAN_EDGE, linewidths=0.8, zorder=4)
    ax.scatter(candidate[x], candidate[y], s=68, c=PURPLE, edgecolors="white", linewidths=0.7, zorder=5)
    if len(reported_low):
        ax.scatter(reported_low[x], reported_low[y], s=38, marker="x", c=ORANGE, linewidths=1.2, zorder=7)
    if len(strict):
        ax.scatter(strict[x], strict[y], s=118, facecolors="none", edgecolors=DARK_PURPLE, linewidths=1.8, zorder=6)
    for _, row in candidate.iterrows():
        ax.annotate(row["formula"], (row[x], row[y]), xytext=(4, 5), textcoords="offset points", fontsize=8, color=DARK_PURPLE)
    ax.axvline(ANALOGUE_THRESHOLD, color="#777777", ls="--", lw=0.8)
    ax.axhline(ANALOGUE_THRESHOLD, color="#777777", ls="--", lw=0.8)
    ax.fill_between([ANALOGUE_THRESHOLD, 1.0], ANALOGUE_THRESHOLD, 1.0, color=PURPLE, alpha=0.06)
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(0.0, 1.01)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("structure-lattice similarity to the same high-zT seed")
    ax.set_ylabel("electronic similarity to that seed")
    ax.set_title(f"{carrier}-type | {title}")
    ax.grid(color="#eeeeee", lw=0.6, zorder=0)


def plot_similarity(frames: dict[str, pd.DataFrame], summaries: dict[str, dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 12.0))
    for column, carrier in enumerate(["n", "p"]):
        frame = frames[carrier]
        scatter_similarity(axes[0, column], frame, "pure", carrier, "pure electronic structure")
        scatter_similarity(axes[1, column], frame, "rich", carrier, "band + fixed-condition transport fingerprint")
        s = summaries[carrier]
        axes[1, column].text(
            0.02,
            0.02,
            f"LOFO AUC: structure {s['lofo_auc_structure_strict']:.2f} | "
            f"band {s['lofo_auc_pure_band']:.2f} | rich electronic {s['lofo_auc_rich_electronic_strict']:.2f}",
            transform=axes[1, column].transAxes,
            fontsize=8.2,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#bbbbbb", alpha=0.9),
        )
    handles = [
        Line2D([], [], marker="o", ls="", color=GREY, label="other material"),
        Line2D([], [], marker="*", ls="", markersize=13, markerfacecolor=CYAN, markeredgecolor=CYAN_EDGE, label="formula-matched experimental zT >= 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor=PURPLE, markeredgecolor="white", label="PF-kL dual screening hit, not high-zT"),
        Line2D([], [], marker="x", ls="", color=ORANGE, label="external formula-level zT report < 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=DARK_PURPLE, markeredgewidth=1.6, label="dual candidate also close in both views"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.012))
    fig.suptitle(
        "Similarity to the same known high-zT analogue: no fitted view weight, no target-driven layout\n"
        "upper-right = close in structure-lattice AND electronic descriptors",
        fontsize=14,
    )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.13, top=0.88, wspace=0.24, hspace=0.31)
    fig.savefig(SIMILARITY_FIGURE, dpi=240, bbox_inches="tight")
    fig.savefig(SIMILARITY_PDF, bbox_inches="tight")
    plt.close(fig)


def plot_and_map(frames: dict[str, pd.DataFrame], summaries: dict[str, dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.3))
    for ax, carrier in zip(axes, ["n", "p"]):
        frame = frames[carrier]
        grey = frame[~frame["known_high_zt"] & ~frame["dual_screening_hit"]]
        seed = frame[frame["known_high_zt"]]
        candidate = frame[frame["dual_screening_hit"]]
        reported_low = frame[frame["reported_low_zt_hit"]]
        strict = frame[frame["strict_analogue"]]
        ax.scatter(grey["and_map_1"], grey["and_map_2"], s=20, c=GREY, alpha=0.55, edgecolors="none", zorder=1)
        ax.scatter(seed["and_map_1"], seed["and_map_2"], s=145, marker="*", facecolors=CYAN, edgecolors=CYAN_EDGE, linewidths=0.8, zorder=5)
        ax.scatter(candidate["and_map_1"], candidate["and_map_2"], s=78, c=PURPLE, edgecolors="white", linewidths=0.8, zorder=6)
        if len(reported_low):
            ax.scatter(reported_low["and_map_1"], reported_low["and_map_2"], s=42, marker="x", c=ORANGE, linewidths=1.25, zorder=8)
        if len(strict):
            ax.scatter(strict["and_map_1"], strict["and_map_2"], s=135, facecolors="none", edgecolors=DARK_PURPLE, linewidths=1.8, zorder=7)
        for _, row in candidate.iterrows():
            seed_row = frame.iloc[int(row["rich_selected_seed_index"])]
            ax.plot(
                [row["and_map_1"], seed_row["and_map_1"]],
                [row["and_map_2"], seed_row["and_map_2"]],
                color=PURPLE,
                alpha=0.38,
                lw=1.0,
                zorder=2,
            )
            ax.annotate(row["formula"], (row["and_map_1"], row["and_map_2"]), xytext=(4, 5), textcoords="offset points", fontsize=8.3, color=DARK_PURPLE)
        for _, row in seed.sort_values("external_zt_max", ascending=False).head(6).iterrows():
            ax.annotate(row["formula"], (row["and_map_1"], row["and_map_2"]), xytext=(4, -10), textcoords="offset points", fontsize=7.7, color=CYAN_EDGE)
        ax.set_xlabel("strict-AND map coordinate 1 (dimensionless)")
        ax.set_ylabel("strict-AND map coordinate 2 (dimensionless)")
        ax.set_title(
            f"{carrier}-type | max(structure rank, electronic rank)\n"
            f"pairwise-rank preservation = {summaries[carrier]['strict_and_map_pairwise_rank_preservation']:.2f}"
        )
        ax.grid(color="#eeeeee", lw=0.6, zorder=0)
    handles = [
        Line2D([], [], marker="o", ls="", color=GREY, label="other material"),
        Line2D([], [], marker="*", ls="", markersize=13, markerfacecolor=CYAN, markeredgecolor=CYAN_EDGE, label="known high-zT formula"),
        Line2D([], [], marker="o", ls="", markerfacecolor=PURPLE, markeredgecolor="white", label="PF-kL dual screening hit"),
        Line2D([], [], marker="x", ls="", color=ORANGE, label="reported formula-level zT < 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=DARK_PURPLE, markeredgewidth=1.6, label="strict analogue"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle(
        "Non-compensatory structure/electronic intersection map\n"
        "the pair distance is the worse neighbour-rank across all descriptor blocks",
        fontsize=14,
    )
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.18, top=0.76, wspace=0.17)
    fig.savefig(AND_FIGURE, dpi=240, bbox_inches="tight")
    fig.savefig(AND_PDF, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    views = pd.read_parquet(VIEWS_PATH).reset_index(drop=True)
    soap = np.load(SOAP_PATH)
    raw = json.load(open(RAW_PATH))
    raw_by_jid = {str(record["jid"]): record for record in raw}
    dual = pd.read_csv(DUAL_PATH)
    zt_by_formula = external_zt_by_formula()

    frames: dict[str, pd.DataFrame] = {}
    summaries: dict[str, dict] = {}
    for carrier in ["n", "p"]:
        frame, summary, _ = build_carrier_frame(
            carrier, views, soap, raw_by_jid, zt_by_formula, dual
        )
        frames[carrier] = frame
        summaries[carrier] = summary

    points = pd.concat(frames.values(), ignore_index=True)
    points.to_csv(POINTS_OUT, index=False)
    candidate_columns = [
        "carrier", "jid", "formula", "PF_jarvis", "kL_exp_300K",
        "PF_percentile", "low_kL_percentile", "dual_desirability",
        "top20_intersection", "pareto", "external_zt_max",
        "pure_structure_similarity", "pure_electronic_similarity", "pure_and_similarity",
        "rich_structure_similarity", "rich_electronic_similarity", "rich_and_similarity",
        "selected_high_zt_seed", "strict_analogue", "external_validation_status",
        "geometry_similarity_to_seed", "composition_lattice_similarity_to_seed",
        "band_similarity_to_seed", "transport_similarity_to_seed",
    ]
    candidates = points[points["dual_screening_hit"]].copy()
    candidates[candidate_columns].sort_values(
        ["carrier", "strict_analogue", "rich_and_similarity", "dual_desirability"],
        ascending=[True, False, False, False],
    ).to_csv(CANDIDATES_OUT, index=False)

    summary_document = {
        "method": {
            "structure": "max(row-rank geometry-only SOAP distance, row-rank composition+lattice-physics distance)",
            "pure_electronic": "row-rank distance in Eg/effective-mass/dielectric/spillage descriptors",
            "rich_electronic": "max(pure-electronic rank, fixed-condition S/sigma/kappa_e fingerprint rank)",
            "joint": "max(structure rank, rich-electronic rank), symmetrized; no fitted block weight",
            "targets_excluded_from_coordinates": ["PF", "experimental kappa_L", "external zT"],
            "important_caveat": "S and sigma algebraically determine PF, so the rich transport fingerprint is analogue visualization, not an independent PF validation.",
        },
        "carrier_results": summaries,
    }
    SUMMARY_OUT.write_text(json.dumps(summary_document, indent=2, ensure_ascii=False), encoding="utf-8")
    plot_similarity(frames, summaries)
    plot_and_map(frames, summaries)

    print(json.dumps(summary_document, indent=2, ensure_ascii=False))
    print("\nNon-high-zT PF-kL screening hits:")
    print(candidates[candidate_columns].to_string(index=False))
    print(f"\nSaved {SIMILARITY_FIGURE}")
    print(f"Saved {AND_FIGURE}")


if __name__ == "__main__":
    main()
