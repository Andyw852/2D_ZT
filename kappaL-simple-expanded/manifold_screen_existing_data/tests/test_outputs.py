import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_target_free_manifold_and_expected_outputs():
    summary = json.loads((ROOT / "outputs/joint_manifold_summary.json").read_text())
    assert summary["no_new_first_principles_or_transport_calculation"] is True
    assert summary["manifold_construction_is_target_free"] is True
    assert summary["n_complete_case_materials"] > 2000
    assert summary["n_seed_formulas"] >= 8
    features = summary["structure_features"] + summary["electronic_features"]
    forbidden = {"PF", "power_factor", "kL", "zt", "zT"}
    assert not any(any(token in name for token in forbidden) for name in features)


def test_candidates_are_unknown_and_figure_exists():
    candidates = pd.read_csv(ROOT / "outputs/manifold_candidate_ranking.csv")
    points = pd.read_csv(ROOT / "outputs/joint_manifold_points.csv")
    assert len(candidates) == 100
    assert candidates["manifold_dual_screen_score"].is_monotonic_decreasing
    known = set(points.loc[points["seed_formula"], "row_id"])
    assert not set(candidates["row_id"]) & known
    figure = ROOT / "figures/joint_structure_electronic_manifold_screen.png"
    assert figure.exists() and figure.stat().st_size > 100_000
