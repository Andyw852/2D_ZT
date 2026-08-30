import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "outputs" / "expanded_simple_space_summary.json"
MEMBERSHIP = ROOT / "outputs" / "expanded_simple_space_membership.csv"
FIGURE = ROOT / "figures" / "te_reference_dual_space_intersection_expanded.png"


def test_expanded_outputs_are_complete_and_consistent():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    membership = pd.read_csv(MEMBERSHIP)

    assert summary["n_expanded_complete_case"] > 8000
    assert len(membership) == summary["n_expanded_complete_case"]
    assert membership["row_id"].is_unique
    assert int(membership["structure_like_top5"].sum()) == summary["n_structure_like"]
    assert int(membership["electronic_like_top5"].sum()) == summary["n_electronic_like"]
    assert int(membership["intersection_top5"].sum()) == summary["n_intersection"]
    assert summary["n_intersection"] > summary["random_expected_intersection"]
    assert summary["n_benchmark_formula_families"] >= 10
    assert FIGURE.stat().st_size > 100_000


def test_intersection_is_the_boolean_overlap():
    membership = pd.read_csv(MEMBERSHIP)
    expected = membership["structure_like_top5"] & membership["electronic_like_top5"]
    assert (membership["intersection_top5"] == expected).all()
    assert (
        membership.loc[membership["robust_intersection"], "intersection_loo_stability"]
        >= 0.80
    ).all()
