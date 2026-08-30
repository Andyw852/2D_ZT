from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_points() -> pd.DataFrame:
    return pd.read_csv(ROOT / "outputs/jarvis_2d_all_points.csv")


def test_all_aligned_2d_materials_are_present_in_each_carrier_overlay() -> None:
    points = load_points()
    assert points.groupby("carrier").size().to_dict() == {"n": 1103, "p": 1103}
    assert points.groupby("carrier")["jid"].nunique().to_dict() == {"n": 1103, "p": 1103}
    assert points.groupby("carrier")["PF_mean"].count().to_dict() == {"n": 806, "p": 803}


def test_global_coordinates_are_reference_free_and_shared() -> None:
    points = load_points()
    summary = json.loads((ROOT / "outputs/jarvis_2d_screen_summary.json").read_text())
    assert summary["coordinates_reference_free"] is True
    assert summary["labels_used_in_coordinates"] == []
    coordinates = ["structure_x", "structure_y", "electronic_x", "electronic_y", "joint_x", "joint_y"]
    assert points[coordinates].notna().all().all()
    n = points[points["carrier"] == "n"].set_index("jid").sort_index()
    p = points[points["carrier"] == "p"].set_index("jid").sort_index()
    assert np.allclose(n[coordinates], p[coordinates])


def test_broad_and_strict_purple_counts_and_thresholds() -> None:
    points = load_points()
    assert points.groupby("carrier")["broad_purple"].sum().to_dict() == {"n": 36, "p": 38}
    assert points.groupby("carrier")["strict_purple"].sum().to_dict() == {"n": 10, "p": 18}
    broad = points[points["broad_purple"]]
    strict = points[points["strict_purple"]]
    assert (broad["PF_percentile"] >= 0.8).all()
    assert (broad["low_kL_surrogate_percentile"] >= 0.8).all()
    assert (strict["PF_percentile"] >= 0.9).all()
    assert (strict["low_kL_surrogate_percentile"] >= 0.9).all()
    assert strict["broad_purple"].all()


def test_candidate_export_matches_broad_screen() -> None:
    candidates = pd.read_csv(ROOT / "outputs/jarvis_2d_purple_candidates.csv")
    assert len(candidates) == 74
    assert candidates["broad_purple"].all()


def test_kappa_surrogate_is_grouped_and_has_ranking_signal() -> None:
    summary = json.loads((ROOT / "outputs/jarvis_2d_screen_summary.json").read_text())
    kappa = summary["kappa_surrogate"]
    assert kappa["training_rows"] == 137
    assert kappa["grouped_cv"].startswith("5-fold by chemical system")
    assert kappa["oof_spearman"] > 0.6
    assert "transferred to JARVIS-2D" in kappa["domain_warning"]


def test_map_local_neighbour_preservation_is_above_random() -> None:
    summary = json.loads((ROOT / "outputs/jarvis_2d_screen_summary.json").read_text())
    preservation = summary["global_maps"]["map_30nn_preservation"]
    assert min(preservation.values()) > 0.5

