import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_strict_and_has_no_view_weight_and_sensitivity():
    summary = json.loads((ROOT / "outputs/strict_and_summary.json").read_text())
    assert summary["view_weight_coefficients"] is None
    assert summary["joint_definition"].startswith("r_AND")
    assert summary["k_sensitivity"] == [15, 30, 50]
    assert summary["manifold_is_target_free"] is True


def test_visual_highlights_are_independent_of_manifold_score():
    points = pd.read_csv(ROOT / "outputs/strict_and_manifold_points.csv")
    assert points["independent_dual_candidate"].sum() == 30
    candidates = pd.read_csv(ROOT / "outputs/strict_and_candidate_ranking.csv")
    assert len(candidates) == 100
    figure = ROOT / "figures/strict_and_joint_manifold.png"
    assert figure.exists() and figure.stat().st_size > 100_000
