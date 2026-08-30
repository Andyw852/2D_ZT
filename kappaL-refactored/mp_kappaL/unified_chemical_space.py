"""Build an asymmetric electronic-aware structural chemical space.

The map is a visualization layer, not a new physical latent-variable claim.
It combines two *superblocks* on the semiconductor complete-case cohort:

1. structural chemistry: composition + composition-blind SOAP, reduced to 30 PCs;
2. recoverable electronics: cross-fitted ExtraTrees predictions of Eg, me*, mh*, eps.

Each superblock is centered and scaled to unit total inertia before fusion.  Thus
splitting the structural input into many columns does not grant it extra block
weight.  The electronic block uses OOF predictions, not in-sample fitted values.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from pymatgen.core import Composition
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config

META_PATH = config.PROC_DIR / "jarvis_structure_electronic_targets.parquet"
COMP_PATH = config.PROC_DIR / "jarvis_composition_magpie_style.npy"
SOAP_PATH = config.PROC_DIR / "jarvis_geometry_soap.npy"
OOF_PATH = config.PROC_DIR / "asymmetric_mapping_oof.parquet"
RAW_PATH = config.EXTERNAL_DATA_DIR / "jarvis_kl" / "jdft_3d-8-18-2021.json"

COORD_OUT = config.PROC_DIR / "unified_chemical_space_coordinates.csv"
BENCHMARK_OUT = config.PROC_DIR / "unified_chemical_space_te_benchmarks.csv"
SENS_OUT = config.PROC_DIR / "unified_chemical_space_weight_sensitivity.csv"
META_OUT = config.PROC_DIR / "unified_chemical_space_metadata.json"
FIG_OUT = config.FIG_DIR / "unified_structure_electronic_space.png"

TARGETS = {
    "positive_gap_log1p_eV": "oof_gap_log1p",
    "electron_mass_log10_me": "oof_log_me",
    "hole_mass_log10_me": "oof_log_mh",
    "dielectric_log10_geomean": "oof_log_eps",
}
STRUCTURE_PCS = 30
K_NEIGHBORS = 15
WEIGHT_RATIOS = (0.5, 1.0, 2.0)

# Reference landmarks, not labels selected or ranked by this model.  Only
# formulas present in the complete-case JARVIS semiconductor cohort are drawn.
TE_BENCHMARK_FORMULAS = (
    "Bi2Te3", "Bi2Se3", "PbTe", "PbSe", "PbS", "SnSe", "GeTe",
    "SnTe", "Cu2Se", "Mg2Si", "SiGe", "TiNiSn", "ZrNiSn", "SrTiO3",
)
BENCHMARK_LABEL_OFFSETS = {
    "Bi2Te3": (-19, 13),
    "Bi2Se3": (7, -13),
    "PbTe": (-18, -13),
    "PbSe": (9, 13),
    "PbS": (2, 14),
    "SnSe": (15, -11),
    "GeTe": (-12, -13),
    "SnTe": (-18, 11),
    "Cu2Se": (8, 9),
    "Mg2Si": (-20, -12),
    "SiGe": (8, 11),
    "TiNiSn": (-5, 12),
    "ZrNiSn": (4, 11),
    "SrTiO3": (-21, 11),
}


def _unit_inertia(block: np.ndarray) -> np.ndarray:
    block = block - block.mean(axis=0, keepdims=True)
    norm = float(np.linalg.norm(block, ord="fro"))
    if norm == 0:
        raise ValueError("Cannot normalize a zero-inertia block")
    return block / norm


def _extract_oof() -> pd.DataFrame:
    oof = pd.read_parquet(OOF_PATH)
    oof = oof[
        (oof["dataset"] == "JARVIS-DFT")
        & (oof["kind"] == "regression")
        & (oof["model"] == "ExtraTrees")
        & (oof["feature_set"] == "C+G")
        & oof["target"].isin(TARGETS)
    ][["target", "row_id", "y_pred"]]
    # Five legacy JIDs are duplicated in the raw 2021 snapshot.  None has valid
    # effective masses and therefore none can enter the complete-case map, but
    # collapse them here so pivoting remains deterministic.
    oof = oof.groupby(["target", "row_id"], as_index=False)["y_pred"].mean()
    wide = oof.pivot(index="row_id", columns="target", values="y_pred")
    wide = wide.rename(columns=TARGETS).dropna(subset=list(TARGETS.values()))
    return wide.reset_index()


def _structural_probes(row_ids: set[str]) -> pd.DataFrame:
    records = json.load(open(RAW_PATH, encoding="utf-8"))
    rows = []
    seen = set()
    for record in records:
        jid = str(record.get("jid"))
        if jid not in row_ids or jid in seen:
            continue
        seen.add(jid)
        atoms = record["atoms"]
        lattice = np.asarray(atoms["lattice_mat"], float)
        lengths = np.linalg.norm(lattice, axis=1)
        angles = np.asarray(atoms.get("angles", (90.0, 90.0, 90.0)), float)
        n_atoms = len(atoms["elements"])
        rows.append({
            "row_id": jid,
            "volume_per_atom": abs(float(np.linalg.det(lattice))) / n_atoms,
            "cell_anisotropy": float(lengths.max() / lengths.min()),
            "angle_distortion": float(np.mean(np.abs(angles - 90.0))),
            "n_atoms_cell": n_atoms,
        })
    result = pd.DataFrame(rows)
    if len(result) != len(row_ids):
        raise ValueError(f"Missing structural probes: {len(result)} != {len(row_ids)}")
    return result


def _mean_neighbour_jaccard(reference: np.ndarray, candidate: np.ndarray) -> float:
    ref = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1).fit(reference)
    cand = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1).fit(candidate)
    a = ref.kneighbors(return_distance=False)[:, 1:]
    b = cand.kneighbors(return_distance=False)[:, 1:]
    values = []
    for left, right in zip(a, b):
        left_set, right_set = set(left), set(right)
        values.append(len(left_set & right_set) / len(left_set | right_set))
    return float(np.mean(values))


def _orient(coords: np.ndarray, gap: np.ndarray, volume: np.ndarray) -> np.ndarray:
    coords = coords.copy()
    if np.corrcoef(coords[:, 0], gap)[0, 1] < 0:
        coords[:, 0] *= -1
    if np.corrcoef(coords[:, 1], volume)[0, 1] < 0:
        coords[:, 1] *= -1
    return coords


def _reduced_formula(formula: str) -> str:
    return Composition(str(formula)).reduced_formula


def _formula_label(formula: str) -> str:
    import re
    return re.sub(r"(\d+)", r"$_{\1}$", formula)


def _benchmark_rows(frame: pd.DataFrame) -> pd.DataFrame:
    lookup = {_reduced_formula(formula): formula for formula in TE_BENCHMARK_FORMULAS}
    tagged = frame.copy()
    tagged["benchmark_formula"] = tagged["formula"].map(_reduced_formula).map(lookup)
    tagged = tagged[tagged["benchmark_formula"].notna()].copy()
    if tagged.empty:
        return tagged

    group = tagged.groupby("benchmark_formula", sort=False)
    tagged["n_polymorphs"] = group["row_id"].transform("size")
    center_x = group["joint_axis_1"].transform("mean")
    center_y = group["joint_axis_2"].transform("mean")
    tagged["distance_to_formula_centroid"] = (
        (tagged["joint_axis_1"] - center_x) ** 2
        + (tagged["joint_axis_2"] - center_y) ** 2
    )
    representative = tagged.groupby("benchmark_formula")[
        "distance_to_formula_centroid"
    ].idxmin()
    tagged["is_label_representative"] = False
    tagged.loc[representative, "is_label_representative"] = True
    return tagged


def _correlation_vectors(coords: np.ndarray, frame: pd.DataFrame, fields: list[str]) -> dict[str, tuple[float, float]]:
    vectors = {}
    for field in fields:
        value = frame[field].to_numpy(float)
        vectors[field] = (
            float(np.corrcoef(coords[:, 0], value)[0, 1]),
            float(np.corrcoef(coords[:, 1], value)[0, 1]),
        )
    return vectors


def _draw_vectors(
    ax: plt.Axes,
    coords: np.ndarray,
    vectors: dict[str, tuple[float, float]],
    labels: dict[str, str],
    color: str,
    label_offsets: dict[str, tuple[float, float]],
) -> None:
    x_radius = np.diff(np.percentile(coords[:, 0], (1, 99)))[0] * 0.32
    y_radius = np.diff(np.percentile(coords[:, 1], (1, 99)))[0] * 0.32
    for field, (vx, vy) in vectors.items():
        ex, ey = vx * x_radius, vy * y_radius
        ax.annotate(
            "", xy=(ex, ey), xytext=(0, 0),
            arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.5, "alpha": .9},
            zorder=5,
        )
        dx, dy = label_offsets.get(field, (0.02, 0.02))
        offset_x = dx * x_radius
        offset_y = dy * y_radius
        ax.text(
            ex + offset_x,
            ey + offset_y,
            labels[field],
            color=color,
            fontsize=9,
            ha="left" if offset_x >= 0 else "right",
            va="bottom" if offset_y >= 0 else "top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .68, "pad": 1.2},
            zorder=6,
        )


def _draw_benchmark_inset(
    ax: plt.Axes,
    frame: pd.DataFrame,
    benchmarks: pd.DataFrame,
) -> None:
    if benchmarks.empty:
        return
    inset = ax.inset_axes([.035, .055, .64, .45])
    bx = benchmarks["joint_axis_1"].to_numpy(float)
    by = benchmarks["joint_axis_2"].to_numpy(float)
    xpad = max(np.ptp(bx) * .09, .002)
    ypad = max(np.ptp(by) * .20, .0015)
    xlim = (bx.min() - xpad, bx.max() + xpad)
    ylim = (by.min() - ypad, by.max() + ypad)
    nearby = frame[
        frame["joint_axis_1"].between(*xlim)
        & frame["joint_axis_2"].between(*ylim)
    ]
    inset.scatter(
        nearby["joint_axis_1"], nearby["joint_axis_2"],
        s=4, c="#a8a8a8", alpha=.24, edgecolors="none", rasterized=True,
    )
    inset.scatter(
        bx, by, marker="*", s=64, c="#00d8e8", edgecolors="#111111",
        linewidths=.55, zorder=4,
    )
    representatives = benchmarks[benchmarks["is_label_representative"]]
    for row in representatives.itertuples(index=False):
        offset = BENCHMARK_LABEL_OFFSETS.get(row.benchmark_formula, (5, 5))
        inset.annotate(
            _formula_label(row.benchmark_formula),
            xy=(row.joint_axis_1, row.joint_axis_2),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.7,
            ha="left" if offset[0] >= 0 else "right",
            va="bottom" if offset[1] >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": .45},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .82, "pad": .35},
            zorder=5,
        )
    inset.set_xlim(xlim)
    inset.set_ylim(ylim)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_title("TE benchmark landmarks (zoom)", fontsize=8.7, pad=2)
    for spine in inset.spines.values():
        spine.set_color("#666666")
        spine.set_linewidth(.65)


def _plot(frame: pd.DataFrame, benchmarks: pd.DataFrame, metadata: dict) -> None:
    coords = frame[["joint_axis_1", "joint_axis_2"]].to_numpy(float)
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 7.7), sharex=True, sharey=True, constrained_layout=True)
    limits = [
        np.percentile(coords[:, 0], (.5, 99.5)),
        np.percentile(coords[:, 1], (.5, 99.5)),
    ]
    pads = [0.13 * np.diff(limit)[0] for limit in limits]

    structure_fields = ["log10_volume_per_atom", "geometry_pc1"]
    structure_labels = {
        "log10_volume_per_atom": "V/atom",
        "geometry_pc1": "SOAP geometry",
    }
    structure_offsets = {
        "log10_volume_per_atom": (-0.03, 0.06),
        "geometry_pc1": (-0.05, 0.13),
    }
    electronic_fields = ["gap_ev", "log10_mass_geomean", "log10_epsilon"]
    electronic_labels = {
        "gap_ev": "Eg",
        "log10_mass_geomean": "m* (e/h)",
        "log10_epsilon": "epsilon",
    }
    electronic_offsets = {
        "gap_ev": (0.03, -0.09),
        "log10_mass_geomean": (0.04, 0.10),
        "log10_epsilon": (-0.04, 0.07),
    }

    color_specs = [
        ("log10_volume_per_atom", "viridis", "log10 volume / atom (Å³)", "Structural lens", structure_fields, structure_labels, structure_offsets, "#c44e20"),
        ("gap_ev", "plasma", "DFT band gap Eg (eV)", "Electronic lens", electronic_fields, electronic_labels, electronic_offsets, "#187b7b"),
    ]
    for ax, (field, cmap, color_label, title, vector_fields, vector_labels, offsets, arrow_color) in zip(axes, color_specs):
        values = frame[field].to_numpy(float)
        lo, hi = np.percentile(values, (2, 98))
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1], c=np.clip(values, lo, hi),
            cmap=cmap, vmin=lo, vmax=hi, s=10, alpha=.58,
            edgecolors="none", rasterized=True,
        )
        ax.axhline(0, color="#dddddd", lw=.7, zorder=0)
        ax.axvline(0, color="#dddddd", lw=.7, zorder=0)
        vectors = _correlation_vectors(coords, frame, vector_fields)
        _draw_vectors(ax, coords, vectors, vector_labels, arrow_color, offsets)
        if not benchmarks.empty:
            ax.scatter(
                benchmarks["joint_axis_1"], benchmarks["joint_axis_2"],
                marker="*", s=92, c="#00d8e8", edgecolors="#111111",
                linewidths=.7, zorder=8,
            )
            ax.legend(
                handles=[Line2D(
                    [0], [0], marker="*", linestyle="none", markersize=10,
                    markerfacecolor="#00d8e8", markeredgecolor="#111111",
                    label="literature TE benchmark",
                )],
                loc="lower right", fontsize=7.5, framealpha=.86,
            )
            _draw_benchmark_inset(ax, frame, benchmarks)
        cb = fig.colorbar(scatter, ax=ax, fraction=.046, pad=.025)
        cb.set_label(color_label)
        ax.set_title(title)
        ax.set_xlabel("joint chemical axis 1")
        ax.set_ylabel("joint chemical axis 2")
        ax.set_xlim(limits[0][0] - pads[0], limits[0][1] + pads[0])
        ax.set_ylim(limits[1][0] - pads[1], limits[1][1] + pads[1])

    fig.suptitle(
        "Unified structure–electronic chemical space\n"
        "Same coordinates: structural lens (left), electronic lens (right) | "
        f"N={len(frame):,}; equal structure/electronic block inertia; 2D variance={metadata['two_dimensional_variance']:.1%}\n"
        "Cyan stars are literature TE landmarks, not model-ranked candidates",
        fontsize=13.2,
    )
    fig.savefig(FIG_OUT, dpi=210, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for path in (META_PATH, COMP_PATH, SOAP_PATH, OOF_PATH, RAW_PATH):
        if not path.exists():
            raise FileNotFoundError(path)
    meta = pd.read_parquet(META_PATH).reset_index(drop=True)
    composition = np.load(COMP_PATH)
    geometry = np.load(SOAP_PATH)
    if len(meta) != len(composition) or len(meta) != len(geometry):
        raise ValueError("JARVIS metadata/features are not row-aligned")

    electronic = _extract_oof()
    core_meta = meta[
        (meta["gap_ev"] > .01)
        & (meta["m_electron"] > 0)
        & (meta["m_hole"] > 0)
        & (meta["epsilon_geo"] > 0)
    ].reset_index(names="source_row")
    if not core_meta["row_id"].is_unique:
        raise ValueError("Complete-case JARVIS JIDs are not unique")
    core = core_meta.merge(
        electronic, on="row_id", how="inner", validate="one_to_one"
    )
    core = core.reset_index(drop=True)
    rows = core["source_row"].to_numpy(int)

    full_structure = np.column_stack([composition[rows], geometry[rows]])
    structure_scaled = StandardScaler().fit_transform(full_structure)
    structure_pca = PCA(n_components=STRUCTURE_PCS, random_state=config.SEED)
    structure_scores = structure_pca.fit_transform(structure_scaled)
    electronic_columns = list(TARGETS.values())
    electronic_scores = StandardScaler().fit_transform(core[electronic_columns])
    structure_block = _unit_inertia(structure_scores)
    electronic_block = _unit_inertia(electronic_scores)

    base_fused = np.column_stack([structure_block, electronic_block])
    joint_pca = PCA(n_components=2, random_state=config.SEED)
    coords = joint_pca.fit_transform(base_fused)

    probes = _structural_probes(set(core["row_id"]))
    core = core.merge(probes, on="row_id", how="left", validate="one_to_one")
    geometry_pc1 = PCA(n_components=1, random_state=config.SEED).fit_transform(
        StandardScaler().fit_transform(geometry[rows])
    )[:, 0]
    coords = _orient(coords, core["gap_ev"].to_numpy(float), core["volume_per_atom"].to_numpy(float))
    core["joint_axis_1"] = coords[:, 0]
    core["joint_axis_2"] = coords[:, 1]
    core["geometry_pc1"] = geometry_pc1
    core["log10_volume_per_atom"] = np.log10(core["volume_per_atom"])
    core["log10_m_electron"] = np.log10(core["m_electron"])
    core["log10_m_hole"] = np.log10(core["m_hole"])
    core["log10_mass_geomean"] = .5 * (
        core["log10_m_electron"] + core["log10_m_hole"]
    )
    core["log10_epsilon"] = np.log10(core["epsilon_geo"])

    base_neighbours = base_fused
    sensitivity_rows = []
    for ratio in WEIGHT_RATIOS:
        fused = np.column_stack([ratio * structure_block, electronic_block])
        pca = PCA(n_components=2, random_state=config.SEED)
        xy = pca.fit_transform(fused)
        sensitivity_rows.append({
            "structure_to_electronic_weight": ratio,
            "structure_inertia_fraction": ratio ** 2 / (ratio ** 2 + 1),
            "two_dimensional_variance": float(pca.explained_variance_ratio_.sum()),
            "knn15_jaccard_vs_equal_weight": _mean_neighbour_jaccard(base_neighbours, fused),
        })
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(SENS_OUT, index=False)

    output_columns = [
        "row_id", "formula", "chemical_system", "joint_axis_1", "joint_axis_2",
        "gap_ev", "m_electron", "m_hole", "epsilon_geo",
        *electronic_columns,
        "volume_per_atom", "cell_anisotropy", "angle_distortion", "n_atoms_cell",
        "geometry_pc1",
    ]
    core[output_columns].to_csv(COORD_OUT, index=False)
    benchmarks = _benchmark_rows(core)
    benchmark_columns = [
        "benchmark_formula", "row_id", "formula", "chemical_system",
        "joint_axis_1", "joint_axis_2", "gap_ev", "m_electron", "m_hole",
        "epsilon_geo", "volume_per_atom", "n_polymorphs",
        "is_label_representative",
    ]
    benchmarks[benchmark_columns].to_csv(BENCHMARK_OUT, index=False)
    metadata = {
        "n_complete_case_semiconductors": len(core),
        "structure_superblock": "162 composition + 147 composition-blind SOAP -> 30 PCs",
        "electronic_superblock": "four cross-fitted ExtraTrees C+G predictions",
        "fusion": "each superblock centered and normalized to unit Frobenius norm",
        "structure_to_electronic_weight": 1.0,
        "structure_inertia_fraction": 0.5,
        "electronic_inertia_fraction": 0.5,
        "structure_pca_variance": float(structure_pca.explained_variance_ratio_.sum()),
        "two_dimensional_variance": float(joint_pca.explained_variance_ratio_.sum()),
        "interpretation": "visualization of recoverable electronic gradients; not a candidate score",
        "te_benchmark_formulas_present": sorted(benchmarks["benchmark_formula"].unique().tolist()),
        "te_benchmark_entries_present": len(benchmarks),
        "te_benchmark_role": "external literature landmarks; not selected or ranked by the model",
    }
    json.dump(metadata, open(META_OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    _plot(core, benchmarks, metadata)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))
    print(sensitivity.to_string(index=False))
    print(f"saved {COORD_OUT}")
    print(f"saved {BENCHMARK_OUT}")
    print(f"saved {SENS_OUT}")
    print(f"saved {FIG_OUT}")


if __name__ == "__main__":
    main()
