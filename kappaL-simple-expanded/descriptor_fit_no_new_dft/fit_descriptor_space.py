"""Fit literature-motivated descriptor blocks using existing local data only.

No DFT, BTE, phonon, or new transport calculation is performed.  The lattice
channel is fitted to the 137 existing starrydata2/JARVIS formula matches with
experimental kL near 300 K.  The electronic channel is fitted to the existing
JARVIS fixed-condition power factors (600 K, 1e20 cm^-3, CRTA).  All electronic
scores are out-of-fold predictions grouped by chemical system.  Lattice scores
for labelled JIDs are also replaced by grouped out-of-fold predictions.
"""
from __future__ import annotations

import importlib.util
import json
import math
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pymatgen.core import Composition
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
WORKSPACE = PARENT.parent
RAW_PATH = (
    WORKSPACE
    / "jarvis_2d_te_atlas"
    / "data"
    / "raw"
    / "external"
    / "jarvis_kl"
    / "jdft_3d-8-18-2021.json"
)
EXPANDED_FEATURES = PARENT / "outputs" / "expanded_simple_space_features.csv"
KL_LABELS = WORKSPACE / "jarvis_2d_te_atlas" / "features" / "kl_verify" / "kl_views.parquet"
ZT_LABELS = (
    WORKSPACE
    / "jarvis_2d_te_atlas"
    / "data"
    / "raw"
    / "external"
    / "starrydata2"
    / "starrydata_zt_sane.csv"
)

OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"
COMPARISON_OUT = OUTPUT_DIR / "descriptor_model_comparison.csv"
SPACE_OUT = OUTPUT_DIR / "cross_validated_descriptor_space.csv"
SUMMARY_OUT = OUTPUT_DIR / "descriptor_space_summary.json"
FIG_OUT = FIGURE_DIR / "cross_validated_structure_electronic_space.png"
FIG_PDF_OUT = FIGURE_DIR / "cross_validated_structure_electronic_space.pdf"

SEED = 20260829
TOP_FRACTION = 0.05
N_SPLITS = 5

COLORS = {
    "other": "#c6c8cc",
    "structure only": "#2f6fed",
    "electronic only": "#f28e2b",
    "intersection": "#8e44ad",
    "external": "#00d8e8",
}

warnings.filterwarnings("ignore", message="No Pauling electronegativity for .*", category=UserWarning)


def _load_parent_module():
    path = PARENT / "run_expanded_dual_space.py"
    spec = importlib.util.spec_from_file_location("expanded_simple_space", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PARENT_MODULE = _load_parent_module()

S0 = (
    "n_elements",
    "composition_entropy",
    "mean_atomic_number",
    "std_atomic_number",
    "log_mean_atomic_mass",
    "mass_coefficient_variation",
    "mean_electronegativity",
    "std_electronegativity",
)
S1 = S0 + (
    "log_density",
    "log_volume_per_atom",
    "log_n_atoms",
    "log_cell_anisotropy",
    "angle_distortion",
)
S2 = S1 + (
    "log1p_bulk_modulus",
    "log1p_shear_modulus",
    "poisson_ratio",
)

E0 = (
    "log1p_gap_ev",
    "log_dielectric_geo",
    "log_dielectric_anisotropy",
)
E1 = E0 + (
    "log1p_mbj_gap",
    "log1p_electron_mass",
    "log1p_hole_mass",
    "log_mass_ratio",
)
E2 = E1 + (
    "log1p_electron_mass_spectral_ratio",
    "log1p_hole_mass_spectral_ratio",
    "log1p_electron_mass_complexity_proxy",
    "log1p_hole_mass_complexity_proxy",
)

BLOCKS = {
    "S0 composition": S0,
    "S1 + geometry": S1,
    "S2 + elastic": S2,
    "E0 gap + dielectric": E0,
    "E1 + effective mass": E1,
    "E2 + mass anisotropy": E2,
}


def finite_positive(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) and result > 0 else float("nan")


def finite_number(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def canon(formula: str) -> str | None:
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return None


def chemical_system(formula: str) -> str:
    try:
        return Composition(str(formula)).chemical_system
    except Exception:
        return "unknown"


def formula_label(formula: str) -> str:
    return re.sub(r"(\d+)", r"$_{\1}$", str(formula))


def add_raw_descriptors(frame: pd.DataFrame, raw_by_jid: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for jid in frame["row_id"]:
        record = raw_by_jid[jid]
        bulk = finite_positive(record.get("bulk_modulus_kv"))
        shear = finite_positive(record.get("shear_modulus_gv"))
        poisson = finite_number(record.get("poisson"))
        mbj = finite_positive(record.get("mbj_bandgap"))
        me = finite_positive(record.get("avg_elec_mass"))
        mh = finite_positive(record.get("avg_hole_mass"))

        masses = record.get("effective_masses_300K") or {}

        def mass_shape(carrier: str) -> tuple[float, float]:
            raw_values = masses.get(carrier, [])
            if not isinstance(raw_values, (list, tuple)):
                return np.nan, np.nan
            values = pd.to_numeric(pd.Series(raw_values), errors="coerce").to_numpy(float)
            values = np.abs(values[np.isfinite(values) & (np.abs(values) > 0)])
            if values.size < 2:
                return np.nan, np.nan
            spectral_ratio = float(values.max() / values.min())
            dos_mass = float(np.prod(values) ** (1.0 / values.size))
            conductivity_mass = float(values.size / np.sum(1.0 / values))
            complexity_proxy = (dos_mass / conductivity_mass) ** 1.5
            return np.log1p(spectral_ratio), np.log1p(complexity_proxy)

        electron_spectral, electron_complexity = mass_shape("n")
        hole_spectral, hole_complexity = mass_shape("p")
        rows.append(
            {
                "row_id": jid,
                "log1p_bulk_modulus": np.log1p(bulk) if np.isfinite(bulk) else np.nan,
                "log1p_shear_modulus": np.log1p(shear) if np.isfinite(shear) else np.nan,
                "poisson_ratio": poisson if 0 < poisson < 0.6 else np.nan,
                "log1p_mbj_gap": np.log1p(mbj) if np.isfinite(mbj) else np.nan,
                "log1p_electron_mass": np.log1p(me) if np.isfinite(me) else np.nan,
                "log1p_hole_mass": np.log1p(mh) if np.isfinite(mh) else np.nan,
                "log_mass_ratio": abs(np.log(me / mh)) if np.isfinite(me) and np.isfinite(mh) else np.nan,
                "log1p_electron_mass_spectral_ratio": electron_spectral,
                "log1p_hole_mass_spectral_ratio": hole_spectral,
                "log1p_electron_mass_complexity_proxy": electron_complexity,
                "log1p_hole_mass_complexity_proxy": hole_complexity,
            }
        )
    return frame.merge(pd.DataFrame(rows), on="row_id", how="left", validate="one_to_one")


def build_expanded_frame(raw_by_jid: dict[str, dict]) -> pd.DataFrame:
    frame = pd.read_csv(EXPANDED_FEATURES)
    frame["canon"] = frame["formula"].map(canon)
    frame["group"] = frame["formula"].map(chemical_system)
    frame = add_raw_descriptors(frame, raw_by_jid)
    frame["electronic_target"] = np.log10(
        1.0 + frame[["power_factor_n_raw", "power_factor_p_raw"]].max(axis=1)
    )
    frame["preferred_carrier"] = np.where(
        frame["power_factor_n_raw"] >= frame["power_factor_p_raw"], "n", "p"
    )
    return frame


def build_lattice_training(raw_by_jid: dict[str, dict]) -> pd.DataFrame:
    labels = pd.read_parquet(KL_LABELS)
    rows = []
    for row in labels.itertuples(index=False):
        record = raw_by_jid.get(str(row.jid))
        if record is None:
            continue
        try:
            values = {
                "row_id": str(row.jid),
                "formula": str(row.formula),
                "canon": canon(str(row.formula)),
                "group": chemical_system(str(row.formula)),
                "kL_300": float(row.kL_300),
            }
            values.update(PARENT_MODULE._composition_features(record["formula"]))
            values.update(PARENT_MODULE._structure_features(record))
            rows.append(values)
        except Exception:
            continue
    frame = pd.DataFrame(rows)
    frame = add_raw_descriptors(frame, raw_by_jid)
    frame["lattice_target"] = -np.log10(frame["kL_300"])
    return frame


def model_factories() -> dict[str, object]:
    elastic = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        RobustScaler(quantile_range=(10.0, 90.0)),
        ElasticNet(alpha=0.03, l1_ratio=0.35, max_iter=20000, random_state=SEED),
    )
    trees = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        ExtraTreesRegressor(
            n_estimators=260,
            min_samples_leaf=3,
            max_features=0.85,
            random_state=SEED,
            n_jobs=-1,
        ),
    )
    return {"ElasticNet": elastic, "ExtraTrees": trees}


def grouped_oof(
    frame: pd.DataFrame,
    features: tuple[str, ...],
    target: str,
    estimator: object,
) -> np.ndarray:
    x = frame.loc[:, features]
    y = frame[target].to_numpy(float)
    groups = frame["group"].astype(str).to_numpy()
    splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    prediction = np.full(len(frame), np.nan)
    for train, test in splitter.split(x, y, groups):
        model = clone(estimator)
        model.fit(x.iloc[train], y[train])
        prediction[test] = model.predict(x.iloc[test])
    if not np.isfinite(prediction).all():
        raise ValueError("Non-finite OOF predictions")
    return prediction


def evaluate_blocks(lattice: pd.DataFrame, electronic: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = []
    predictions: dict[tuple[str, str], np.ndarray] = {}
    for block, features in BLOCKS.items():
        channel = "structure" if block.startswith("S") else "electronic"
        frame = lattice if channel == "structure" else electronic
        target = "lattice_target" if channel == "structure" else "electronic_target"
        for model_name, estimator in model_factories().items():
            pred = grouped_oof(frame, features, target, estimator)
            predictions[(block, model_name)] = pred
            rho = float(spearmanr(frame[target], pred).statistic)
            rows.append(
                {
                    "channel": channel,
                    "block": block,
                    "model": model_name,
                    "n_samples": len(frame),
                    "n_features": len(features),
                    "spearman_oof": rho,
                    "r2_oof": float(r2_score(frame[target], pred)),
                    "mae_oof": float(mean_absolute_error(frame[target], pred)),
                }
            )
    result = pd.DataFrame(rows)
    selected = {}
    for channel in ("structure", "electronic"):
        subset = result[result["channel"] == channel]
        row = subset.sort_values(["spearman_oof", "r2_oof"], ascending=False).iloc[0]
        selected[channel] = {
            "block": str(row.block),
            "model": str(row.model),
            "features": list(BLOCKS[str(row.block)]),
            "spearman_oof": float(row.spearman_oof),
            "r2_oof": float(row.r2_oof),
            "oof_prediction": predictions[(str(row.block), str(row.model))],
        }
    return result, selected


def fit_selected_and_build_space(
    lattice: pd.DataFrame,
    electronic: pd.DataFrame,
    selected: dict,
) -> pd.DataFrame:
    structure_spec = selected["structure"]
    electronic_spec = selected["electronic"]
    factories = model_factories()

    structure_model = clone(factories[structure_spec["model"]])
    structure_model.fit(lattice[structure_spec["features"]], lattice["lattice_target"])
    structure_prediction = structure_model.predict(electronic[structure_spec["features"]])
    labelled_oof = dict(zip(lattice["row_id"], structure_spec["oof_prediction"]))
    structure_prediction = np.asarray(
        [labelled_oof.get(jid, value) for jid, value in zip(electronic["row_id"], structure_prediction)],
        dtype=float,
    )
    electronic_prediction = np.asarray(electronic_spec["oof_prediction"], dtype=float)

    out = electronic.copy()
    out["predicted_low_kL_score_raw"] = structure_prediction
    out["predicted_electronic_score_raw"] = electronic_prediction
    out["structure_score_percentile"] = pd.Series(structure_prediction).rank(pct=True).to_numpy()
    out["electronic_score_percentile"] = pd.Series(electronic_prediction).rank(pct=True).to_numpy()
    out["dual_score"] = np.sqrt(
        out["structure_score_percentile"] * out["electronic_score_percentile"]
    )
    out["structure_top5"] = out["structure_score_percentile"] >= 1.0 - TOP_FRACTION
    out["electronic_top5"] = out["electronic_score_percentile"] >= 1.0 - TOP_FRACTION
    out["intersection_top5"] = out["structure_top5"] & out["electronic_top5"]
    category = np.full(len(out), "other", dtype=object)
    category[out["structure_top5"] & ~out["electronic_top5"]] = "structure only"
    category[~out["structure_top5"] & out["electronic_top5"]] = "electronic only"
    category[out["intersection_top5"]] = "intersection"
    out["category"] = category
    return out


def add_external_zt_labels(frame: pd.DataFrame) -> pd.DataFrame:
    zt = pd.read_csv(ZT_LABELS)
    zt = zt[(zt["T_max"] >= 200.0) & (zt["T_max"] <= 1500.0)].copy()
    zt["canon"] = zt["composition"].map(canon)
    zt = zt.dropna(subset=["canon", "zt_max"])
    by_formula = zt.groupby("canon", as_index=False).agg(
        external_zt_max=("zt_max", "max"),
        external_zt_temperature_max=("T_max", "max"),
        external_zt_samples=("sample_id", "nunique"),
    )
    out = frame.merge(by_formula, on="canon", how="left", validate="many_to_one")
    out["external_high_zt_formula"] = out["external_zt_max"].fillna(-np.inf) >= 1.0
    return out


def external_formula_recall(frame: pd.DataFrame, threshold: float) -> tuple[int, int, float]:
    high = frame[frame["external_high_zt_formula"]]
    formulas = high["canon"].dropna().unique()
    recovered = high[
        (high["structure_score_percentile"] >= threshold)
        & (high["electronic_score_percentile"] >= threshold)
    ]["canon"].dropna().unique()
    return len(formulas), len(recovered), len(recovered) / len(formulas) if len(formulas) else float("nan")


def plot(comparison: pd.DataFrame, frame: pd.DataFrame, summary: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18.2, 6.5), constrained_layout=True)

    ax = axes[0]
    best_by_block = (
        comparison.sort_values("spearman_oof", ascending=False)
        .drop_duplicates("block")
        .set_index("block")
        .loc[list(BLOCKS)]
        .reset_index()
    )
    y = np.arange(len(best_by_block))[::-1]
    colors = ["#2f6fed" if channel == "structure" else "#f28e2b" for channel in best_by_block["channel"]]
    ax.barh(y, best_by_block["spearman_oof"], color=colors, alpha=0.88)
    for yy, row in zip(y, best_by_block.itertuples(index=False)):
        ax.text(
            max(row.spearman_oof, 0.0) + 0.012,
            yy,
            f"ρ={row.spearman_oof:.2f}  {row.model}",
            va="center",
            fontsize=8.5,
        )
    ax.axvline(0, color="#777777", lw=0.7)
    ax.set_yticks(y, best_by_block["block"])
    ax.set_xlim(min(-0.1, best_by_block["spearman_oof"].min() - 0.05), 0.82)
    ax.set_xlabel("grouped-CV Spearman")
    ax.set_title("Descriptor-block comparison\nbest fixed model per block")

    order = ["other", "structure only", "electronic only", "intersection"]
    sizes = {"other": 5, "structure only": 22, "electronic only": 22, "intersection": 52}
    alpha = {"other": 0.13, "structure only": 0.82, "electronic only": 0.82, "intersection": 0.97}
    external = frame[frame["external_high_zt_formula"]]

    for ax, limits, title in (
        (axes[1], (0.0, 1.005), "CV-calibrated descriptor score space\nall 9,029 materials"),
        (axes[2], (0.70, 1.005), "Upper-right decision region\ntop 30% shown"),
    ):
        for category in order:
            subset = frame[frame["category"] == category]
            ax.scatter(
                subset["electronic_score_percentile"],
                subset["structure_score_percentile"],
                s=sizes[category],
                color=COLORS[category],
                alpha=alpha[category],
                edgecolors="#222222" if category == "intersection" else "none",
                linewidths=0.45 if category == "intersection" else 0,
                rasterized=category == "other",
                zorder=1 if category == "other" else 3,
            )
        ax.scatter(
            external["electronic_score_percentile"],
            external["structure_score_percentile"],
            marker="*",
            s=86,
            color=COLORS["external"],
            edgecolors="#111111",
            linewidths=0.65,
            zorder=7,
        )
        ax.axvline(0.95, color=COLORS["electronic only"], lw=1.0, ls="--")
        ax.axhline(0.95, color=COLORS["structure only"], lw=1.0, ls="--")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_xlabel("electronic PF prediction percentile (OOF)")
        ax.set_ylabel("predicted low-kL score percentile")
        ax.set_title(title)

    labels = (
        external.sort_values(["external_zt_max", "dual_score"], ascending=False)
        .drop_duplicates("canon")
        .head(10)
    )
    for row in labels.itertuples(index=False):
        if row.electronic_score_percentile < 0.70 or row.structure_score_percentile < 0.70:
            continue
        axes[2].annotate(
            formula_label(row.canon),
            (row.electronic_score_percentile, row.structure_score_percentile),
            xytext=(-5, 5),
            textcoords="offset points",
            fontsize=7.8,
            ha="right",
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "#555555"},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.3},
            zorder=8,
        )

    legend = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["other"], markeredgecolor="none", markersize=6, label="other"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["structure only"], markeredgecolor="none", markersize=7, label="low-kL only"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["electronic only"], markeredgecolor="none", markersize=7, label="high-PF only"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLORS["intersection"], markeredgecolor="#222222", markersize=8, label=f"top-5% intersection (N={summary['n_intersection_top5']})"),
        Line2D([0], [0], marker="*", linestyle="none", markerfacecolor=COLORS["external"], markeredgecolor="#111111", markersize=10, label="starrydata2 formula match, max zT≥1"),
    ]
    axes[1].legend(handles=legend, loc="lower left", fontsize=8.0, framealpha=0.88)
    axes[2].text(
        0.705,
        0.708,
        f"high-zT formula recall in top 5% × 5%: {summary['external_high_zt_top5_recovered_formulas']}/"
        f"{summary['external_high_zt_unique_formulas']}\n"
        f"top 10% × 10%: {summary['external_high_zt_top10_recovered_formulas']}/"
        f"{summary['external_high_zt_unique_formulas']}",
        fontsize=8.8,
        va="bottom",
    )
    fig.suptitle(
        "Literature-motivated descriptor fitting with existing data only\n"
        "x is grouped-OOF; y is CV-selected from 137 experimental labels (labelled JIDs use OOF); zT stars are formula-level checks",
        fontsize=14,
    )
    fig.savefig(FIG_OUT, dpi=230, bbox_inches="tight")
    fig.savefig(FIG_PDF_OUT, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    raw = json.load(open(RAW_PATH, encoding="utf-8"))
    raw_by_jid = {str(record["jid"]): record for record in raw}

    electronic = build_expanded_frame(raw_by_jid)
    lattice = build_lattice_training(raw_by_jid)
    comparison, selected = evaluate_blocks(lattice, electronic)
    space = fit_selected_and_build_space(lattice, electronic, selected)
    space = add_external_zt_labels(space)

    n5, recovered5, recall5 = external_formula_recall(space, 0.95)
    n10, recovered10, recall10 = external_formula_recall(space, 0.90)
    summary = {
        "no_new_first_principles_calculation": True,
        "n_materials_in_space": int(len(space)),
        "n_lattice_training_labels": int(len(lattice)),
        "lattice_label": "starrydata2 experimental kL near 300 K, formula-matched to JARVIS",
        "electronic_label": "max(n,p) JARVIS power factor at 600 K and 1e20 cm^-3, CRTA",
        "selected_structure_block": {k: v for k, v in selected["structure"].items() if k != "oof_prediction"},
        "selected_electronic_block": {k: v for k, v in selected["electronic"].items() if k != "oof_prediction"},
        "n_structure_top5": int(space["structure_top5"].sum()),
        "n_electronic_top5": int(space["electronic_top5"].sum()),
        "n_intersection_top5": int(space["intersection_top5"].sum()),
        "external_high_zt_unique_formulas": int(n5),
        "external_high_zt_matching_jids": int(space["external_high_zt_formula"].sum()),
        "external_high_zt_top5_recovered_formulas": int(recovered5),
        "external_high_zt_top5_formula_recall": float(recall5),
        "external_high_zt_top10_recovered_formulas": int(recovered10),
        "external_high_zt_top10_formula_recall": float(recall10),
        "external_validation_warning": (
            "starrydata2 zT is matched by reduced formula; polymorph, doping, dimensionality, "
            "microstructure and measurement conditions may differ from the JARVIS structure"
        ),
    }

    comparison.to_csv(COMPARISON_OUT, index=False)
    selected_columns = [
        "row_id",
        "formula",
        "canon",
        "chemical_system",
        "preferred_carrier",
        "power_factor_n_raw",
        "power_factor_p_raw",
        "predicted_low_kL_score_raw",
        "predicted_electronic_score_raw",
        "structure_score_percentile",
        "electronic_score_percentile",
        "dual_score",
        "structure_top5",
        "electronic_top5",
        "intersection_top5",
        "category",
        "external_zt_max",
        "external_zt_temperature_max",
        "external_zt_samples",
        "external_high_zt_formula",
    ]
    space[selected_columns].to_csv(SPACE_OUT, index=False)
    with open(SUMMARY_OUT, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    plot(comparison, space, summary)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nmodel comparison")
    print(comparison.sort_values(["channel", "spearman_oof"], ascending=[True, False]).to_string(index=False))
    print("\nintersection")
    print(
        space.loc[space["intersection_top5"], ["row_id", "formula", "dual_score", "external_zt_max"]]
        .sort_values("dual_score", ascending=False)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
