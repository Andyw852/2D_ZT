"""Reference-free JARVIS-2D structure/electronic manifold and broad screen.

All 1,103 materials with structure and electronic descriptors determine the
three global maps.  Existing n/p PF labels and a cross-domain experimental-kL
structure surrogate are overlaid only after the coordinates are fixed.

No DFT, BTE, phonon, or new transport calculation is performed.
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
from scipy.linalg import orthogonal_procrustes
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, pairwise_distances, r2_score
from sklearn.metrics.pairwise import cosine_distances
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import RobustScaler
from umap import UMAP


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[1]
ATLAS = WORKSPACE / "jarvis_2d_te_atlas"

SOAP_PATH = ATLAS / "features/structure/geometry_soap_v1.parquet"
MAGPIE_PATH = ATLAS / "features/structure/composition_magpie.parquet"
ELECTRONIC_PATH = ATLAS / "features/electronic/electronic_features_v1.parquet"
TRANSPORT_PATH = ATLAS / "features/transport/{carrier}_transport_tensor_features.parquet"
RAW_2D_PATH = ATLAS / "data/raw/jarvis/dft_2d_snapshot.json"

KL_VIEWS_PATH = ATLAS / "features/kl_verify/kl_views.parquet"
KL_SOAP_PATH = ATLAS / "data/processed/kl_soap_geo.npy"
RAW_3D_PATH = ATLAS / "data/raw/external/jarvis_kl/jdft_3d-8-18-2021.json"
ZT_PATH = ATLAS / "data/raw/external/starrydata2/starrydata_zt_sane.csv"

OUTPUT_DIR = HERE / "outputs"
FIGURE_DIR = HERE / "figures"
POINTS_OUT = OUTPUT_DIR / "jarvis_2d_all_points.csv"
CANDIDATES_OUT = OUTPUT_DIR / "jarvis_2d_purple_candidates.csv"
SUMMARY_OUT = OUTPUT_DIR / "jarvis_2d_screen_summary.json"
FIGURE_OUT = FIGURE_DIR / "jarvis_2d_global_structure_electronic_screen.png"
PDF_OUT = FIGURE_DIR / "jarvis_2d_global_structure_electronic_screen.pdf"

SEED = 20260830
BROAD_QUANTILE = 0.80
STRICT_QUANTILE = 0.90

GREY = "#c9cdd2"
CYAN = "#00cbd8"
CYAN_EDGE = "#008c96"
PURPLE = "#9b59b6"
DARK_PURPLE = "#5b247a"
ORANGE = "#e07a16"


def canon(formula: object) -> str | None:
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return None


def chemical_system(formula: object) -> str:
    try:
        return "-".join(sorted(str(element) for element in Composition(str(formula)).elements))
    except Exception:
        return str(formula)


def finite(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return np.nan
    return result if np.isfinite(result) and result > -99998 else np.nan


def load_2d_identity() -> pd.DataFrame:
    records = json.load(open(RAW_2D_PATH))
    rows = []
    for record in records:
        attributes = record.get("attributes", record)
        jid = attributes.get("_jarvis_jid", attributes.get("jid"))
        formula = attributes.get("_jarvis_formula", attributes.get("formula"))
        if jid:
            rows.append({"jid": str(jid), "formula": formula, "canon": canon(formula)})
    return pd.DataFrame(rows).drop_duplicates("jid")


def external_zt() -> pd.DataFrame:
    zt = pd.read_csv(ZT_PATH)
    zt = zt[(zt["T_max"] >= 200.0) & (zt["T_max"] <= 1500.0)].copy()
    zt["canon"] = zt["composition"].map(canon)
    zt = zt.dropna(subset=["canon", "zt_max"])
    return zt.groupby("canon", as_index=False).agg(
        external_zt_max=("zt_max", "max"),
        external_zt_samples=("sample_id", "nunique"),
    )


PROPERTY_NAMES = [
    "electronegativity", "atomic_mass", "atomic_radius", "row", "group",
    "ionization_energy", "electron_affinity", "Z",
]


def element_property(symbol: str, name: str) -> float:
    element = Element(symbol)
    value = {
        "electronegativity": element.X,
        "atomic_mass": element.atomic_mass,
        "atomic_radius": element.atomic_radius,
        "row": element.row,
        "group": element.group,
        "ionization_energy": element.ionization_energy,
        "electron_affinity": element.electron_affinity,
        "Z": element.Z,
    }[name]
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def magpie_from_species(species: list[str]) -> dict[str, float]:
    features: dict[str, float] = {}
    for name in PROPERTY_NAMES:
        values = np.asarray([element_property(symbol, name) for symbol in species], dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            for statistic in ["mean", "std", "min", "max", "range"]:
                features[f"{name}_{statistic}"] = np.nan
            continue
        features[f"{name}_mean"] = float(values.mean())
        features[f"{name}_std"] = float(values.std())
        features[f"{name}_min"] = float(values.min())
        features[f"{name}_max"] = float(values.max())
        features[f"{name}_range"] = float(values.max() - values.min())
    return features


def l2_normalize(values: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.where(norm > 0, norm, 1.0)


def row_rank(distance: np.ndarray) -> np.ndarray:
    n = len(distance)
    order = np.argsort(distance, axis=1, kind="mergesort")
    rank = np.empty((n, n), dtype=np.float32)
    rank[np.arange(n)[:, None], order] = np.arange(n, dtype=np.float32)[None, :]
    rank /= max(1, n - 1)
    np.fill_diagonal(rank, 0.0)
    return rank


def numeric_distance(frame: pd.DataFrame) -> np.ndarray:
    values = SimpleImputer(strategy="median", add_indicator=True).fit_transform(frame.to_numpy(float))
    values = RobustScaler(quantile_range=(10.0, 90.0)).fit_transform(values)
    values = np.clip(values, -8.0, 8.0)
    return pairwise_distances(values, metric="euclidean") / math.sqrt(max(1, values.shape[1]))


def symmetric_rank_distance(distance: np.ndarray) -> np.ndarray:
    rank = row_rank(distance)
    result = np.maximum(rank, rank.T)
    np.fill_diagonal(result, 0.0)
    return result


def normalized_map(distance: np.ndarray) -> np.ndarray:
    coordinates = UMAP(
        n_components=2,
        n_neighbors=35,
        min_dist=0.18,
        metric="precomputed",
        random_state=SEED,
        init="random",
    ).fit_transform(distance)
    coordinates -= coordinates.mean(axis=0, keepdims=True)
    scale = np.sqrt(np.mean(np.sum(coordinates**2, axis=1)))
    return coordinates / max(scale, 1e-12)


def align(reference: np.ndarray, moving: np.ndarray) -> np.ndarray:
    rotation, _ = orthogonal_procrustes(moving, reference)
    return moving @ rotation


def sampled_rank_preservation(distance: np.ndarray, coordinates: np.ndarray) -> float:
    rng = np.random.default_rng(SEED)
    n = len(distance)
    i = rng.integers(0, n, size=100_000)
    j = rng.integers(0, n, size=100_000)
    keep = i != j
    embedded = np.linalg.norm(coordinates[i[keep]] - coordinates[j[keep]], axis=1)
    return float(spearmanr(distance[i[keep], j[keep]], embedded).statistic)


def knn_preservation(distance: np.ndarray, coordinates: np.ndarray, k: int = 30) -> float:
    embedded_distance = pairwise_distances(coordinates)
    original_neighbours = np.argsort(distance, axis=1, kind="mergesort")[:, 1 : k + 1]
    embedded_neighbours = np.argsort(embedded_distance, axis=1, kind="mergesort")[:, 1 : k + 1]
    overlaps = [
        len(np.intersect1d(original_neighbours[row], embedded_neighbours[row], assume_unique=False)) / k
        for row in range(len(distance))
    ]
    return float(np.mean(overlaps))


def build_global_maps(master: pd.DataFrame) -> tuple[pd.DataFrame, dict, dict[str, np.ndarray]]:
    soap_columns = {
        cutoff: [column for column in master.columns if column.startswith(f"soap{cutoff}_mean_")]
        for cutoff in [4, 6, 8]
    }
    magpie_columns = [
        column for column in pd.read_parquet(MAGPIE_PATH).columns if column != "jid"
    ]
    electronic_columns = [
        column for column in pd.read_parquet(ELECTRONIC_PATH).columns if column != "jid"
    ]

    structure_subviews = []
    for cutoff in [4, 6, 8]:
        soap = l2_normalize(master[soap_columns[cutoff]].to_numpy(float))
        structure_subviews.append(symmetric_rank_distance(cosine_distances(soap)))
    structure_subviews.append(symmetric_rank_distance(numeric_distance(master[magpie_columns])))
    structure_distance = np.maximum.reduce(structure_subviews)
    electronic_distance = symmetric_rank_distance(numeric_distance(master[electronic_columns]))
    joint_distance = np.maximum(structure_distance, electronic_distance)

    maps = {
        "structure": normalized_map(structure_distance),
        "electronic": normalized_map(electronic_distance),
        "joint": normalized_map(joint_distance),
    }
    maps["electronic"] = align(maps["structure"], maps["electronic"])
    maps["joint"] = align(maps["structure"], maps["joint"])

    out = master.copy()
    for name, coordinates in maps.items():
        out[f"{name}_x"] = coordinates[:, 0]
        out[f"{name}_y"] = coordinates[:, 1]
    summary = {
        "n_global_materials": int(len(out)),
        "structure_blocks": {"SOAP_r4": 147, "SOAP_r6": 147, "SOAP_r8": 147, "Magpie": 40},
        "electronic_features": electronic_columns,
        "map_pairwise_rank_preservation": {
            "structure": sampled_rank_preservation(structure_distance, maps["structure"]),
            "electronic": sampled_rank_preservation(electronic_distance, maps["electronic"]),
            "joint": sampled_rank_preservation(joint_distance, maps["joint"]),
        },
        "map_30nn_preservation": {
            "structure": knn_preservation(structure_distance, maps["structure"]),
            "electronic": knn_preservation(electronic_distance, maps["electronic"]),
            "joint": knn_preservation(joint_distance, maps["joint"]),
        },
    }
    return out, summary, {
        "structure": structure_distance,
        "electronic": electronic_distance,
        "joint": joint_distance,
    }


def train_kl_surrogate(master: pd.DataFrame) -> tuple[np.ndarray, dict]:
    kl = pd.read_parquet(KL_VIEWS_PATH).reset_index(drop=True)
    raw = json.load(open(RAW_3D_PATH))
    raw_by_jid = {str(record["jid"]): record for record in raw}
    magpie_columns = [column for column in pd.read_parquet(MAGPIE_PATH).columns if column != "jid"]

    training_magpie = []
    for jid in kl["jid"].astype(str):
        species = list((raw_by_jid[jid].get("atoms") or {}).get("elements", []))
        training_magpie.append(magpie_from_species(species))
    training_magpie = pd.DataFrame(training_magpie).reindex(columns=magpie_columns)
    target_magpie = master[magpie_columns]

    soap_train = l2_normalize(np.load(KL_SOAP_PATH).astype(float))
    soap_columns = [column for column in master.columns if column.startswith("soap6_mean_")]
    soap_target = l2_normalize(master[soap_columns].to_numpy(float))

    imputer = SimpleImputer(strategy="median", add_indicator=False)
    magpie_train_array = imputer.fit_transform(training_magpie.to_numpy(float))
    magpie_target_array = imputer.transform(target_magpie.to_numpy(float))
    x_train = np.hstack([soap_train, magpie_train_array])
    x_target = np.hstack([soap_target, magpie_target_array])
    y = np.log10(kl["kL_300"].to_numpy(float))
    groups = kl["formula"].map(chemical_system).to_numpy()

    model_parameters = dict(
        n_estimators=600,
        max_features=0.65,
        min_samples_leaf=2,
        random_state=SEED,
        n_jobs=-1,
    )
    oof = np.full(len(y), np.nan)
    splitter = GroupKFold(n_splits=5)
    for train, test in splitter.split(x_train, y, groups):
        model = ExtraTreesRegressor(**model_parameters)
        model.fit(x_train[train], y[train])
        oof[test] = model.predict(x_train[test])
    model = ExtraTreesRegressor(**model_parameters)
    model.fit(x_train, y)
    prediction = np.clip(model.predict(x_target), y.min(), y.max())

    summary = {
        "training_rows": int(len(y)),
        "target_rows": int(len(prediction)),
        "target": "log10 experimental kappa_L near 300 K",
        "grouped_cv": "5-fold by chemical system",
        "oof_r2": float(r2_score(y, oof)),
        "oof_spearman": float(spearmanr(y, oof).statistic),
        "oof_mae_log10": float(mean_absolute_error(y, oof)),
        "domain_warning": "trained on 137 formula-matched JARVIS-3D structures and transferred to JARVIS-2D; ranking surrogate only",
    }
    return prediction, summary


def add_screen(master: pd.DataFrame, carrier: str) -> pd.DataFrame:
    tensor = pd.read_parquet(Path(str(TRANSPORT_PATH).format(carrier=carrier)))
    transport_columns = [
        "jid", "S_mean", "S_std", "S_abs_mean", "S_relative_spread",
        "log_sigma_mean", "sigma_anisotropy_log", "log_kappa_e_mean",
        "kappa_e_anisotropy_log", "PF_mean", "PF_anisotropy_log",
    ]
    out = master.merge(tensor[transport_columns], on="jid", how="left", validate="one_to_one")
    covered = out["PF_mean"].notna() & (out["PF_mean"] > 0)
    out["PF_percentile"] = np.nan
    out.loc[covered, "PF_percentile"] = out.loc[covered, "PF_mean"].rank(pct=True)
    # A lower predicted log-kL is better.  Rank only within the carrier-covered
    # set so the percentile has the same denominator as PF.
    out["low_kL_surrogate_percentile"] = np.nan
    out.loc[covered, "low_kL_surrogate_percentile"] = (
        1.0 - out.loc[covered, "predicted_log10_kL"].rank(pct=True) + 1.0 / int(covered.sum())
    )
    out["broad_purple"] = covered & (
        (out["PF_percentile"] >= BROAD_QUANTILE)
        & (out["low_kL_surrogate_percentile"] >= BROAD_QUANTILE)
    )
    out["strict_purple"] = covered & (
        (out["PF_percentile"] >= STRICT_QUANTILE)
        & (out["low_kL_surrogate_percentile"] >= STRICT_QUANTILE)
    )
    out["known_high_zt"] = out["external_zt_max"].fillna(-np.inf) >= 1.0
    out["known_high_zt_and_broad"] = out["known_high_zt"] & out["broad_purple"]
    out["reported_low_zt_purple"] = out["broad_purple"] & out["external_zt_max"].notna() & ~out["known_high_zt"]
    out["carrier"] = carrier
    return out


def draw_panel(ax, frame: pd.DataFrame, prefix: str, title: str, preservation: float) -> None:
    x, y = f"{prefix}_x", f"{prefix}_y"
    high = frame[frame["known_high_zt"]]
    broad = frame[frame["broad_purple"] & ~frame["known_high_zt"]]
    strict = frame[frame["strict_purple"] & ~frame["known_high_zt"]]
    overlap = frame[frame["known_high_zt_and_broad"]]
    low_report = frame[frame["reported_low_zt_purple"]]

    ax.scatter(frame[x], frame[y], s=10, c=GREY, alpha=0.32, edgecolors="none", rasterized=True, zorder=1)
    ax.scatter(broad[x], broad[y], s=38, c=PURPLE, alpha=0.88, edgecolors="white", linewidths=0.35, zorder=4)
    ax.scatter(strict[x], strict[y], s=70, c=DARK_PURPLE, edgecolors="white", linewidths=0.55, zorder=5)
    if len(low_report):
        ax.scatter(low_report[x], low_report[y], s=28, marker="x", c=ORANGE, linewidths=0.9, zorder=7)
    ax.scatter(high[x], high[y], s=105, marker="*", c=CYAN, edgecolors=CYAN_EDGE, linewidths=0.7, zorder=6)
    if len(overlap):
        ax.scatter(overlap[x], overlap[y], s=150, facecolors="none", edgecolors=DARK_PURPLE, linewidths=1.5, zorder=8)

    if prefix == "joint":
        offsets = [(3, 5), (3, -11), (-34, 5), (-34, -11)]
        labelled = strict.sort_values(
            ["PF_percentile", "low_kL_surrogate_percentile"], ascending=False
        ).head(4)
        for (_, row), offset in zip(labelled.iterrows(), offsets):
            ax.annotate(row["formula"], (row[x], row[y]), xytext=offset, textcoords="offset points", fontsize=6.7, color=DARK_PURPLE)
    for _, row in high.sort_values("external_zt_max", ascending=False).head(3).iterrows():
        ax.annotate(row["formula"], (row[x], row[y]), xytext=(3, -9), textcoords="offset points", fontsize=6.5, color=CYAN_EDGE)

    ax.set_xlabel("global coordinate 1 (dimensionless)")
    ax.set_ylabel("global coordinate 2 (dimensionless)")
    ax.set_title(f"{title}\n30-NN preservation = {preservation:.2f}")
    ax.grid(color="#eeeeee", lw=0.5, zorder=0)


def plot(frames: dict[str, pd.DataFrame], map_summary: dict, kappa_summary: dict) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17.0, 10.6))
    panels = [
        ("structure", "global structure space"),
        ("electronic", "global electronic-structure space"),
        ("joint", "strict-AND common space"),
    ]
    for row, carrier in enumerate(["n", "p"]):
        for column, (prefix, title) in enumerate(panels):
            draw_panel(
                axes[row, column], frames[carrier], prefix,
                f"{carrier}-type | {title}",
                map_summary["map_30nn_preservation"][prefix],
            )

    handles = [
        Line2D([], [], marker="o", ls="", color=GREY, label="all other JARVIS-2D materials"),
        Line2D([], [], marker="o", ls="", markerfacecolor=PURPLE, markeredgecolor="white", label="top-20% PF AND top-20% low-kL surrogate"),
        Line2D([], [], marker="o", ls="", markerfacecolor=DARK_PURPLE, markeredgecolor="white", label="strict top-10% AND top-10%"),
        Line2D([], [], marker="*", ls="", markersize=12, markerfacecolor=CYAN, markeredgecolor=CYAN_EDGE, label="formula-matched experimental zT >= 1"),
        Line2D([], [], marker="x", ls="", color=ORANGE, label="purple with external formula-level zT < 1"),
        Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=DARK_PURPLE, markeredgewidth=1.5, label="known high-zT also passes broad screen"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        "All usable JARVIS-2D data: reference-free structure/electronic common space\n"
        "all 1,103 materials set the coordinates; n/p transport covers 806/803",
        fontsize=15,
    )
    fig.text(
        0.5, 0.082,
        f"Low-kL colour is a cross-domain structure surrogate (137 experimental training rows; grouped-CV Spearman={kappa_summary['oof_spearman']:.2f}), not a 2D experimental value",
        ha="center", fontsize=9.3, color="#444444",
    )
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.17, top=0.86, wspace=0.22, hspace=0.34)
    fig.savefig(FIGURE_OUT, dpi=240, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    identity = load_2d_identity()
    soap = pd.read_parquet(SOAP_PATH)
    magpie = pd.read_parquet(MAGPIE_PATH)
    electronic = pd.read_parquet(ELECTRONIC_PATH)
    master = soap.merge(magpie, on="jid", validate="one_to_one")
    master = master.merge(electronic, on="jid", validate="one_to_one")
    master = master.merge(identity, on="jid", how="left", validate="one_to_one")
    master = master.merge(external_zt(), on="canon", how="left", validate="many_to_one")
    master, map_summary, _ = build_global_maps(master)

    predicted_log_kappa, kappa_summary = train_kl_surrogate(master)
    master["predicted_log10_kL"] = predicted_log_kappa
    frames = {carrier: add_screen(master, carrier) for carrier in ["n", "p"]}

    points = pd.concat(frames.values(), ignore_index=True)
    output_columns = [
        "carrier", "jid", "formula", "canon", "external_zt_max", "known_high_zt",
        "PF_mean", "PF_percentile", "predicted_log10_kL", "low_kL_surrogate_percentile",
        "broad_purple", "strict_purple", "known_high_zt_and_broad", "reported_low_zt_purple",
        "S_mean", "S_std", "S_abs_mean", "S_relative_spread", "log_sigma_mean",
        "sigma_anisotropy_log", "log_kappa_e_mean", "kappa_e_anisotropy_log", "PF_anisotropy_log",
        "structure_x", "structure_y", "electronic_x", "electronic_y", "joint_x", "joint_y",
    ]
    points[output_columns].to_csv(POINTS_OUT, index=False)
    candidates = points[points["broad_purple"]].copy()
    candidates[output_columns].sort_values(
        ["carrier", "strict_purple", "PF_percentile", "low_kL_surrogate_percentile"],
        ascending=[True, False, False, False],
    ).to_csv(CANDIDATES_OUT, index=False)

    carrier_summary = {}
    for carrier, frame in frames.items():
        carrier_summary[carrier] = {
            "transport_covered": int(frame["PF_mean"].notna().sum()),
            "broad_purple_top20_intersection": int(frame["broad_purple"].sum()),
            "strict_purple_top10_intersection": int(frame["strict_purple"].sum()),
            "high_zt_formula_rows": int(frame["known_high_zt"].sum()),
            "high_zt_also_broad": int(frame["known_high_zt_and_broad"].sum()),
            "reported_low_zt_purple": int(frame["reported_low_zt_purple"].sum()),
        }
    summary = {
        "scope": "all 1,103 aligned JARVIS-2D structure+electronic rows; all available n/p transport rows",
        "coordinates_reference_free": True,
        "labels_used_in_coordinates": [],
        "purple_definition": {
            "broad": "PF percentile >= 0.80 AND low-kL surrogate percentile >= 0.80",
            "strict": "PF percentile >= 0.90 AND low-kL surrogate percentile >= 0.90",
        },
        "global_maps": map_summary,
        "kappa_surrogate": kappa_summary,
        "carriers": carrier_summary,
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    plot(frames, map_summary, kappa_summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nPurple candidates by carrier:")
    print(candidates.groupby(["carrier", "strict_purple"]).size().to_string())
    print(f"Saved {FIGURE_OUT}")


if __name__ == "__main__":
    main()
