from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_expected_sample_sizes_and_labels() -> None:
    points = pd.read_csv(ROOT / "outputs/rich_descriptor_points.csv")
    assert points.groupby("carrier").size().to_dict() == {"n": 85, "p": 85}
    assert points.groupby("carrier")["known_high_zt"].sum().to_dict() == {"n": 11, "p": 11}


def test_strict_analogue_is_noncompensatory() -> None:
    candidates = pd.read_csv(ROOT / "outputs/rich_descriptor_candidates.csv")
    strict = candidates[candidates["strict_analogue"]]
    assert (strict["rich_structure_similarity"] >= 0.75).all()
    assert (strict["rich_electronic_similarity"] >= 0.75).all()
    assert (strict["rich_and_similarity"] == strict[[
        "rich_structure_similarity", "rich_electronic_similarity"
    ]].min(axis=1)).all()


def test_targets_are_declared_excluded_from_coordinates() -> None:
    summary = json.loads((ROOT / "outputs/rich_descriptor_summary.json").read_text())
    excluded = set(summary["method"]["targets_excluded_from_coordinates"])
    assert excluded == {"PF", "experimental kappa_L", "external zT"}
    assert "algebraically determine PF" in summary["method"]["important_caveat"]


def test_every_selected_seed_is_high_zt_and_not_self_for_seed_rows() -> None:
    points = pd.read_csv(ROOT / "outputs/rich_descriptor_points.csv")
    high_formulas = set(points.loc[points["known_high_zt"], "canon"])
    assert set(points["selected_high_zt_seed"]).issubset(high_formulas)
    seed_rows = points[points["known_high_zt"]]
    assert (seed_rows["canon"] != seed_rows["selected_high_zt_seed"]).all()


def test_global_common_space_is_reference_free_and_finite() -> None:
    summary = json.loads((ROOT / "outputs/global_common_space_summary.json").read_text())
    points = pd.read_csv(ROOT / "outputs/global_common_space_points.csv")
    assert summary["reference_free"] is True
    assert summary["labels_used_in_coordinates"] == []
    assert points.groupby("carrier").size().to_dict() == {"n": 85, "p": 85}
    coordinate_columns = [
        "structure_x", "structure_y", "pure_electronic_x", "pure_electronic_y",
        "pure_joint_and_x", "pure_joint_and_y",
    ]
    assert points[coordinate_columns].notna().all().all()


def test_global_joint_space_excludes_transport_targets() -> None:
    summary = json.loads((ROOT / "outputs/global_common_space_summary.json").read_text())
    assert set(summary["electronic_space_excludes"]) == {"S", "sigma", "kappa_e", "PF", "zT"}
    assert summary["joint_distance"].startswith("max(")
