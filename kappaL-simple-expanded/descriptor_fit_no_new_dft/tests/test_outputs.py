import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "outputs" / "descriptor_space_summary.json"
COMPARISON = ROOT / "outputs" / "descriptor_model_comparison.csv"
SPACE = ROOT / "outputs" / "cross_validated_descriptor_space.csv"
FIGURE = ROOT / "figures" / "cross_validated_structure_electronic_space.png"


def test_descriptor_space_outputs_are_consistent():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    comparison = pd.read_csv(COMPARISON)
    space = pd.read_csv(SPACE)

    assert summary["no_new_first_principles_calculation"] is True
    assert len(space) == summary["n_materials_in_space"] == 9029
    assert summary["n_lattice_training_labels"] == 137
    assert len(comparison) == 12
    assert space["row_id"].is_unique
    assert int(space["intersection_top5"].sum()) == summary["n_intersection_top5"]
    assert space["structure_score_percentile"].between(0, 1).all()
    assert space["electronic_score_percentile"].between(0, 1).all()
    assert FIGURE.stat().st_size > 100_000


def test_no_target_columns_are_selected_as_descriptors():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    selected = (
        summary["selected_structure_block"]["features"]
        + summary["selected_electronic_block"]["features"]
    )
    forbidden = {"kL_300", "lattice_target", "power_factor_n_raw", "power_factor_p_raw", "electronic_target"}
    assert forbidden.isdisjoint(selected)
