"""Asymmetric structure -> electronic-property decisive experiment.

This module implements the pre-registered comparison described in
``reports/03_asymmetric_structure_to_transport_plan.md``:

* JARVIS records are analysed in their native material-level structures, so
  polymorphs are retained and no MP/JARVIS formula replication is needed.
* A Magpie-style composition baseline (element fractions + weighted elemental
  statistics) is compared with composition + composition-blind geometric SOAP.
* Outer folds hold out complete chemical systems.  Every feature-set comparison
  uses identical folds and out-of-fold predictions.
* The MP Snyder target is kept only as a leakage/model-internal positive control;
  the small material-resolved experimental kappa_L cohort is reported separately.

The output is diagnostic.  It does not calculate absolute weighted mobility or B,
because the local electronic transport data contain sigma/tau rather than sigma.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
from jarvis.core.atoms import Atoms as JarvisAtoms
from pymatgen.core import Composition, Element
from scipy import stats
from sklearn.base import clone
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config

JARVIS_RAW = (
    config.EXTERNAL_DATA_DIR / "jarvis_kl" / "jdft_3d-8-18-2021.json"
)
SEED = config.SEED
N_SPLITS = 5
N_BOOT = 500
GAP_THRESHOLD_EV = 0.01

JARVIS_META = config.PROC_DIR / "jarvis_structure_electronic_targets.parquet"
JARVIS_COMP = config.PROC_DIR / "jarvis_composition_magpie_style.npy"
JARVIS_SOAP = config.PROC_DIR / "jarvis_geometry_soap.npy"
FEATURE_META = config.PROC_DIR / "asymmetric_feature_metadata.json"
SUMMARY_OUT = config.PROC_DIR / "asymmetric_mapping_summary.csv"
FOLDS_OUT = config.PROC_DIR / "asymmetric_mapping_folds.csv"
OOF_OUT = config.PROC_DIR / "asymmetric_mapping_oof.parquet"
FIG_OUT = config.FIG_DIR / "asymmetric_mapping_decisive_test.png"


def _finite(value) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    if not np.isfinite(value) or value <= -99998:
        return np.nan
    return value


def _element_value(el: Element, name: str) -> float:
    value = getattr(el, name, None)
    if value is None:
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


ELEMENT_PROPERTIES = (
    "Z",
    "atomic_mass",
    "row",
    "group",
    "X",
    "atomic_radius",
    "mendeleev_no",
)
STATISTICS = ("mean", "std", "min", "max", "range", "missing_fraction")


def composition_features_from_amounts(
    amounts: dict[str, float],
) -> tuple[np.ndarray, list[str], str]:
    """Element fractions plus weighted elemental-property statistics.

    The construction is deliberately independent of compound DFT labels.  It is
    Magpie-style rather than a dependency on matminer, which is not in the pinned
    environment.
    """
    total = float(sum(amounts.values()))
    fractions = {symbol: float(n) / total for symbol, n in amounts.items()}
    fraction_vector = np.zeros(118, dtype=np.float32)
    for symbol, fraction in fractions.items():
        fraction_vector[Element(symbol).Z - 1] = fraction

    values: list[float] = []
    names = [f"fraction_Z{z}" for z in range(1, 119)]
    for prop in ELEMENT_PROPERTIES:
        prop_values = []
        prop_weights = []
        missing = 0.0
        for symbol, fraction in fractions.items():
            value = _element_value(Element(symbol), prop)
            if np.isfinite(value):
                prop_values.append(value)
                prop_weights.append(fraction)
            else:
                missing += fraction
        if prop_values:
            x = np.asarray(prop_values, float)
            w = np.asarray(prop_weights, float)
            w /= w.sum()
            mean = float(np.sum(w * x))
            std = float(np.sqrt(np.sum(w * (x - mean) ** 2)))
            summary = (mean, std, float(x.min()), float(x.max()),
                       float(x.max() - x.min()), float(missing))
        else:
            summary = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        values.extend(summary)
        names.extend(f"{prop}_{stat}" for stat in STATISTICS)

    f = np.asarray(list(fractions.values()), float)
    entropy = float(-(f * np.log(f)).sum())
    values.extend((float(len(fractions)), entropy))
    names.extend(("n_elements", "composition_entropy"))
    chemical_system = "-".join(sorted(fractions))
    return np.concatenate([fraction_vector, np.asarray(values, np.float32)]), names, chemical_system


def composition_features(elements: list[str]) -> tuple[np.ndarray, list[str], str]:
    return composition_features_from_amounts(dict(Counter(elements)))


def _jarvis_ase(record: dict) -> Atoms:
    atoms = JarvisAtoms.from_dict(record["atoms"])
    return Atoms(
        symbols=["X"] * atoms.num_atoms,
        positions=np.asarray(atoms.cart_coords, float),
        cell=np.asarray(atoms.lattice_mat, float),
        pbc=True,
    )


def build_jarvis_features(rebuild: bool = False) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    cached = all(path.exists() for path in (JARVIS_META, JARVIS_COMP, JARVIS_SOAP))
    if cached and not rebuild:
        return (
            pd.read_parquet(JARVIS_META),
            np.load(JARVIS_COMP),
            np.load(JARVIS_SOAP),
        )

    if not JARVIS_RAW.exists():
        raise FileNotFoundError(f"Missing local JARVIS data: {JARVIS_RAW}")
    records = json.load(open(JARVIS_RAW, encoding="utf-8"))
    rows: list[dict] = []
    comp_rows: list[np.ndarray] = []
    geometry_chunks: list[np.ndarray] = []
    structure_chunk: list[Atoms] = []
    comp_names: list[str] | None = None
    soap = SOAP(
        species=config.SOAP_SPECIES,
        periodic=True,
        r_cut=config.SOAP_R_CUT,
        n_max=config.SOAP_N_MAX,
        l_max=config.SOAP_L_MAX,
        sigma=config.SOAP_SIGMA,
        average=config.SOAP_AVERAGE,
    )
    for record in records:
        elements = list(record["atoms"]["elements"])
        comp, names, chemical_system = composition_features(elements)
        if comp_names is None:
            comp_names = names
        comp_rows.append(comp)
        structure_chunk.append(_jarvis_ase(record))
        eps = [_finite(record.get(key)) for key in ("epsx", "epsy", "epsz")]
        eps_geo = (
            float(np.exp(np.mean(np.log(eps))))
            if all(np.isfinite(eps)) and all(value > 0 for value in eps)
            else np.nan
        )
        rows.append({
            "row_id": record.get("jid"),
            "formula": record.get("formula"),
            "chemical_system": chemical_system,
            "n_atoms": len(elements),
            "gap_ev": _finite(record.get("optb88vdw_bandgap")),
            "m_electron": _finite(record.get("avg_elec_mass")),
            "m_hole": _finite(record.get("avg_hole_mass")),
            "epsilon_geo": eps_geo,
            "reference": record.get("reference"),
        })
        if len(structure_chunk) == 4000:
            geometry_chunks.append(
                np.asarray(soap.create(structure_chunk, n_jobs=-1), dtype=np.float32)
            )
            print(f"SOAP: {len(rows):,}/{len(records):,}")
            structure_chunk = []
    if structure_chunk:
        geometry_chunks.append(
            np.asarray(soap.create(structure_chunk, n_jobs=-1), dtype=np.float32)
        )
    geometry = np.vstack(geometry_chunks)
    composition = np.asarray(comp_rows, dtype=np.float32)
    metadata = pd.DataFrame(rows)
    if len(metadata) != len(composition) or len(metadata) != len(geometry):
        raise RuntimeError("JARVIS feature rows are not aligned")
    metadata.to_parquet(JARVIS_META, index=False)
    np.save(JARVIS_COMP, composition)
    np.save(JARVIS_SOAP, geometry)
    json.dump({
        "composition_features": comp_names,
        "composition_shape": list(composition.shape),
        "geometry_shape": list(geometry.shape),
        "geometry_definition": "composition-blind SOAP on native JARVIS optimized cells",
        "soap": {
            "r_cut": config.SOAP_R_CUT,
            "n_max": config.SOAP_N_MAX,
            "l_max": config.SOAP_L_MAX,
            "sigma": config.SOAP_SIGMA,
            "average": config.SOAP_AVERAGE,
        },
    }, open(FEATURE_META, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    return metadata, composition, geometry


def _regression_models() -> dict[str, object]:
    return {
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=96,
            min_samples_leaf=3,
            max_features=0.55,
            n_jobs=-1,
            random_state=SEED,
        ),
    }


def _classification_models() -> dict[str, object]:
    return {
        "Logistic": make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=1.0,
                class_weight="balanced",
                max_iter=2000,
                random_state=SEED,
            ),
        ),
        "ExtraTrees": ExtraTreesClassifier(
            n_estimators=96,
            min_samples_leaf=3,
            max_features=0.55,
            class_weight="balanced",
            n_jobs=-1,
            random_state=SEED,
        ),
    }


def _score(kind: str, y: np.ndarray, prediction: np.ndarray) -> float:
    if kind == "classification":
        return float(roc_auc_score(y, prediction))
    return float(r2_score(y, prediction))


def _bootstrap_delta(
    kind: str,
    y: np.ndarray,
    pred_c: np.ndarray,
    pred_cg: np.ndarray,
    groups: np.ndarray,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    unique = np.unique(groups)
    by_group = {group: np.flatnonzero(groups == group) for group in unique}
    deltas = []
    for _ in range(N_BOOT):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([by_group[group] for group in sampled])
        if kind == "classification" and np.unique(y[indices]).size < 2:
            continue
        deltas.append(
            _score(kind, y[indices], pred_cg[indices])
            - _score(kind, y[indices], pred_c[indices])
        )
    if not deltas:
        return np.nan, np.nan
    return tuple(float(x) for x in np.quantile(deltas, (0.025, 0.975)))


def evaluate_target(
    *,
    dataset: str,
    target: str,
    kind: str,
    y: np.ndarray,
    groups: np.ndarray,
    row_ids: np.ndarray,
    feature_sets: dict[str, np.ndarray],
    summary_rows: list[dict],
    fold_rows: list[dict],
    oof_rows: list[pd.DataFrame],
) -> None:
    models = _classification_models() if kind == "classification" else _regression_models()
    splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    splits = list(splitter.split(np.zeros(len(y)), y, groups))
    metric = "ROC-AUC" if kind == "classification" else "R2"
    for model_name, base_model in models.items():
        predictions: dict[str, np.ndarray] = {}
        fold_ids = np.full(len(y), -1, dtype=np.int8)
        for feature_name, x in feature_sets.items():
            pred = np.full(len(y), np.nan, dtype=float)
            for fold, (train, test) in enumerate(splits):
                model = clone(base_model)
                model.fit(x[train], y[train])
                if kind == "classification":
                    pred[test] = model.predict_proba(x[test])[:, 1]
                else:
                    pred[test] = model.predict(x[test])
                fold_ids[test] = fold
                fold_score = _score(kind, y[test], pred[test])
                fold_rows.append({
                    "dataset": dataset,
                    "target": target,
                    "kind": kind,
                    "model": model_name,
                    "feature_set": feature_name,
                    "fold": fold,
                    "n_test": len(test),
                    "n_test_groups": len(np.unique(groups[test])),
                    "metric": metric,
                    "score": fold_score,
                    "mae": (
                        np.nan if kind == "classification"
                        else float(mean_absolute_error(y[test], pred[test]))
                    ),
                    "spearman": (
                        np.nan if kind == "classification"
                        else float(stats.spearmanr(y[test], pred[test]).statistic)
                    ),
                })
            predictions[feature_name] = pred

        score_c = _score(kind, y, predictions["C"])
        score_cg = _score(kind, y, predictions["C+G"])
        ci_lo, ci_hi = _bootstrap_delta(
            kind, y, predictions["C"], predictions["C+G"], groups,
            seed=SEED + sum(ord(char) for char in f"{dataset}-{target}-{model_name}"),
        )
        for feature_name, pred in predictions.items():
            score = _score(kind, y, pred)
            summary_rows.append({
                "dataset": dataset,
                "target": target,
                "kind": kind,
                "model": model_name,
                "feature_set": feature_name,
                "n": len(y),
                "n_groups": len(np.unique(groups)),
                "metric": metric,
                "score": score,
                "mae": (
                    np.nan if kind == "classification"
                    else float(mean_absolute_error(y, pred))
                ),
                "spearman": (
                    np.nan if kind == "classification"
                    else float(stats.spearmanr(y, pred).statistic)
                ),
                "delta_CG_minus_C": score_cg - score_c if feature_name == "C+G" else np.nan,
                "delta_ci_lo": ci_lo if feature_name == "C+G" else np.nan,
                "delta_ci_hi": ci_hi if feature_name == "C+G" else np.nan,
                "grouping": "held-out chemical_system",
            })
            oof_rows.append(pd.DataFrame({
                "dataset": dataset,
                "target": target,
                "kind": kind,
                "model": model_name,
                "feature_set": feature_name,
                "row_id": row_ids,
                "chemical_system": groups,
                "fold": fold_ids,
                "y_true": y,
                "y_pred": pred,
            }))


def electronic_experiment(
    meta: pd.DataFrame,
    composition: np.ndarray,
    geometry: np.ndarray,
    summary_rows: list[dict],
    fold_rows: list[dict],
    oof_rows: list[pd.DataFrame],
) -> None:
    x_all = np.column_stack([composition, geometry])
    target_specs = [
        (
            "gap_nonmetal",
            "classification",
            meta["gap_ev"].notna().to_numpy(),
            lambda frame: (frame["gap_ev"].to_numpy(float) > GAP_THRESHOLD_EV).astype(int),
        ),
        (
            "positive_gap_log1p_eV",
            "regression",
            (meta["gap_ev"] > GAP_THRESHOLD_EV).to_numpy(),
            lambda frame: np.log1p(frame["gap_ev"].to_numpy(float)),
        ),
        (
            "electron_mass_log10_me",
            "regression",
            (meta["m_electron"] > 0).to_numpy(),
            lambda frame: np.log10(frame["m_electron"].to_numpy(float)),
        ),
        (
            "hole_mass_log10_me",
            "regression",
            (meta["m_hole"] > 0).to_numpy(),
            lambda frame: np.log10(frame["m_hole"].to_numpy(float)),
        ),
        (
            "dielectric_log10_geomean",
            "regression",
            (meta["epsilon_geo"] > 0).to_numpy(),
            lambda frame: np.log10(frame["epsilon_geo"].to_numpy(float)),
        ),
    ]
    for target, kind, mask, transform in target_specs:
        frame = meta.loc[mask]
        evaluate_target(
            dataset="JARVIS-DFT",
            target=target,
            kind=kind,
            y=transform(frame),
            groups=frame["chemical_system"].to_numpy(str),
            row_ids=frame["row_id"].to_numpy(str),
            feature_sets={"C": composition[mask], "C+G": x_all[mask]},
            summary_rows=summary_rows,
            fold_rows=fold_rows,
            oof_rows=oof_rows,
        )


def _mp_feature_blocks() -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet").reset_index(drop=True)
    row_ids = np.load(config.PROC_DIR / "row_index.npy", allow_pickle=True)
    if not np.array_equal(meta["material_id"].to_numpy(), row_ids):
        raise ValueError("MP metadata and feature arrays are not aligned")
    geometry = np.load(config.PROC_DIR / "soap_geo.npy")
    comp_rows = []
    groups = []
    for formula in meta["formula"]:
        composition = Composition(str(formula))
        features, _, group = composition_features_from_amounts(
            composition.get_el_amt_dict()
        )
        comp_rows.append(features)
        groups.append(group)
    composition = np.asarray(comp_rows, np.float32)
    avg_mass = composition[:, 118 + ELEMENT_PROPERTIES.index("atomic_mass") * len(STATISTICS)]
    physics = np.column_stack([
        np.log1p(meta["bulk_vrh"].to_numpy(float)),
        np.log1p(meta["shear_vrh"].to_numpy(float)),
        np.log1p(meta["debye"].to_numpy(float)),
        np.log1p(meta["density"].to_numpy(float)),
        np.log1p(meta["nsites"].to_numpy(float)),
        np.log1p(meta["v_long"].to_numpy(float)),
        np.log1p(meta["v_trans"].to_numpy(float)),
        np.log1p(avg_mass),
    ]).astype(np.float32)
    meta["chemical_system"] = groups
    return meta, {
        "C": composition,
        "C+G": np.column_stack([composition, geometry]),
        "P_kappa": physics,
    }


def lattice_experiment(
    summary_rows: list[dict],
    fold_rows: list[dict],
    oof_rows: list[pd.DataFrame],
) -> None:
    meta, features = _mp_feature_blocks()
    valid = np.isfinite(meta["snyder_acoustic"]) & (meta["snyder_acoustic"] > 0)
    evaluate_target(
        dataset="MP",
        target="Snyder_300K_model_log10",
        kind="regression",
        y=np.log10(meta.loc[valid, "snyder_acoustic"].to_numpy(float)),
        groups=meta.loc[valid, "chemical_system"].to_numpy(str),
        row_ids=meta.loc[valid, "material_id"].to_numpy(str),
        feature_sets={name: x[valid] for name, x in features.items()},
        summary_rows=summary_rows,
        fold_rows=fold_rows,
        oof_rows=oof_rows,
    )

    targets = pd.read_parquet(config.PROC_DIR / "kappa_L_targets.parquet")
    exp = targets[
        (targets["method"] == "experimental")
        & (targets["match_quality"] == "unique_formula_to_mpid")
        & targets["material_id"].notna()
        & (targets["kappa_L"] > 0)
    ][["material_id", "kappa_L"]].drop_duplicates("material_id")
    joined = meta[["material_id", "chemical_system"]].merge(
        exp, on="material_id", how="inner", validate="one_to_one"
    )
    position = pd.Series(np.arange(len(meta)), index=meta["material_id"])
    idx = position.loc[joined["material_id"]].to_numpy(int)
    evaluate_target(
        dataset="Starrydata2×MP",
        target="experimental_kappaL_log10",
        kind="regression",
        y=np.log10(joined["kappa_L"].to_numpy(float)),
        groups=joined["chemical_system"].to_numpy(str),
        row_ids=joined["material_id"].to_numpy(str),
        feature_sets={name: x[idx] for name, x in features.items()},
        summary_rows=summary_rows,
        fold_rows=fold_rows,
        oof_rows=oof_rows,
    )


TARGET_LABELS = {
    "positive_gap_log1p_eV": "positive gap",
    "electron_mass_log10_me": "electron mass",
    "hole_mass_log10_me": "hole mass",
    "dielectric_log10_geomean": "dielectric",
}


def plot_results(summary: pd.DataFrame) -> None:
    electronic = summary[
        (summary["dataset"] == "JARVIS-DFT")
        & (summary["kind"] == "regression")
        & summary["target"].isin(TARGET_LABELS)
    ].copy()
    order = list(TARGET_LABELS)
    y_pos = np.arange(len(order))[::-1]

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), constrained_layout=True)
    colors = {"C": "#8c8c8c", "C+G": "#157f7a", "P_kappa": "#c66a19"}
    markers = {"Ridge": "o", "ExtraTrees": "s"}

    for ax, model in zip(axes[0], ("Ridge", "ExtraTrees")):
        subset = electronic[electronic["model"] == model]
        for y, target in zip(y_pos, order):
            rows = subset[subset["target"] == target].set_index("feature_set")
            if not {"C", "C+G"}.issubset(rows.index):
                continue
            x0, x1 = rows.loc["C", "score"], rows.loc["C+G", "score"]
            ax.plot([x0, x1], [y, y], color="#bdbdbd", lw=1.5, zorder=1)
            ax.scatter(x0, y, s=55, color=colors["C"], marker="o", label="composition C" if y == y_pos[0] else None)
            ax.scatter(x1, y, s=64, color=colors["C+G"], marker="s", label="C + geometry G" if y == y_pos[0] else None)
        ax.axvline(0, color="#b0b0b0", lw=.8)
        ax.set_yticks(y_pos, [TARGET_LABELS[target] for target in order])
        ax.set_xlabel("out-of-fold R²")
        ax.set_title(f"Electronic targets — {model}")
        ax.grid(axis="x", color="#e6e6e6", lw=.7)
        ax.legend(frameon=False, fontsize=9, loc="lower right")

    ax = axes[1, 0]
    offsets = {"Ridge": -0.10, "ExtraTrees": 0.10}
    for model in ("Ridge", "ExtraTrees"):
        subset = electronic[(electronic["model"] == model) & (electronic["feature_set"] == "C+G")]
        xs, ys, lo, hi = [], [], [], []
        for y, target in zip(y_pos, order):
            row = subset[subset["target"] == target].iloc[0]
            xs.append(row["delta_CG_minus_C"])
            ys.append(y + offsets[model])
            lo.append(row["delta_CG_minus_C"] - row["delta_ci_lo"])
            hi.append(row["delta_ci_hi"] - row["delta_CG_minus_C"])
        ax.errorbar(xs, ys, xerr=[lo, hi], fmt=markers[model], ms=6, capsize=3,
                    lw=1.3, label=model)
    ax.axvline(0, color="#333333", lw=1)
    ax.set_yticks(y_pos, [TARGET_LABELS[target] for target in order])
    ax.set_xlabel("ΔR² = R²(C+G) − R²(C), 95% group bootstrap CI")
    ax.set_title("Does geometry add information beyond composition?")
    ax.grid(axis="x", color="#e6e6e6", lw=.7)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 1]
    lattice = summary[
        summary["target"].isin(("Snyder_300K_model_log10", "experimental_kappaL_log10"))
    ].copy()
    x_positions = {
        ("Snyder_300K_model_log10", "Ridge"): 0,
        ("Snyder_300K_model_log10", "ExtraTrees"): 1,
        ("experimental_kappaL_log10", "Ridge"): 3,
        ("experimental_kappaL_log10", "ExtraTrees"): 4,
    }
    for (target, model), x in x_positions.items():
        rows = lattice[(lattice["target"] == target) & (lattice["model"] == model)].set_index("feature_set")
        vals = [rows.loc[name, "score"] for name in ("C", "C+G", "P_kappa")]
        ax.plot([x] * 3, vals, color="#cfcfcf", lw=1.2, zorder=1)
        for name, value in zip(("C", "C+G", "P_kappa"), vals):
            ax.scatter(x, value, s=62, color=colors[name],
                       marker={"C": "o", "C+G": "s", "P_kappa": "^"}[name],
                       label=name if x == 0 else None, zorder=2)
    ax.axhline(0, color="#888888", lw=.8)
    ax.set_xticks([.5, 3.5], ["Snyder model\nN=12,156", "experimental κL\nN=59"])
    ax.set_ylabel("out-of-fold R²")
    ax.set_title("Lattice channel: leakage control vs low-power experiment")
    ax.grid(axis="y", color="#e6e6e6", lw=.7)
    ax.legend(frameon=False, fontsize=9, ncol=3, loc="lower left")

    classification = summary[
        (summary["target"] == "gap_nonmetal")
        & (summary["feature_set"].isin(("C", "C+G")))
    ]
    notes = []
    for model in ("Logistic", "ExtraTrees"):
        rows = classification[classification["model"] == model].set_index("feature_set")
        if {"C", "C+G"}.issubset(rows.index):
            notes.append(f"{model} nonmetal AUC {rows.loc['C','score']:.3f}→{rows.loc['C+G','score']:.3f}")
    fig.suptitle(
        "Decisive test: structure-to-electronic prediction under held-out chemical systems\n"
        + "; ".join(notes),
        fontsize=14,
    )
    fig.savefig(FIG_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rebuild-features", action="store_true",
        help="Recompute cached JARVIS composition and geometry descriptors.",
    )
    args = parser.parse_args()
    config.PROC_DIR.mkdir(parents=True, exist_ok=True)
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)
    metadata, composition, geometry = build_jarvis_features(args.rebuild_features)
    print(f"JARVIS feature cohort: {len(metadata):,}; C={composition.shape}; G={geometry.shape}")

    summary_rows: list[dict] = []
    fold_rows: list[dict] = []
    oof_rows: list[pd.DataFrame] = []
    electronic_experiment(metadata, composition, geometry, summary_rows, fold_rows, oof_rows)
    lattice_experiment(summary_rows, fold_rows, oof_rows)

    summary = pd.DataFrame(summary_rows)
    folds = pd.DataFrame(fold_rows)
    oof = pd.concat(oof_rows, ignore_index=True)
    summary.to_csv(SUMMARY_OUT, index=False)
    folds.to_csv(FOLDS_OUT, index=False)
    oof.to_parquet(OOF_OUT, index=False)
    plot_results(summary)

    display = summary[summary["feature_set"] == "C+G"][
        ["dataset", "target", "model", "n", "metric", "score",
         "delta_CG_minus_C", "delta_ci_lo", "delta_ci_hi"]
    ]
    print(display.to_string(index=False))
    print(f"saved {SUMMARY_OUT}")
    print(f"saved {FOLDS_OUT}")
    print(f"saved {OOF_OUT}")
    print(f"saved {FIG_OUT}")


if __name__ == "__main__":
    main()
