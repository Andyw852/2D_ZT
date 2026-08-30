import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_consensus_outputs_are_present_and_aligned():
    points = pd.read_csv(ROOT / "outputs" / "consensus_audit_points.csv")
    candidates = pd.read_csv(ROOT / "outputs" / "consensus_candidates.csv")
    summary = json.loads(
        (ROOT / "outputs" / "consensus_audit_summary.json").read_text(
            encoding="utf-8"
        )
    )

    assert len(points) == summary["n_materials"]
    assert points["row_id"].is_unique
    assert int(points["consensus_candidate"].sum()) == summary["n_consensus_candidates"]
    assert len(candidates) == summary["n_consensus_candidates"]
    assert set(candidates["row_id"]) == set(
        points.loc[points["consensus_candidate"], "row_id"]
    )


def test_consensus_candidates_obey_noncompensatory_rule():
    points = pd.read_csv(ROOT / "outputs" / "consensus_audit_points.csv")
    summary = json.loads(
        (ROOT / "outputs" / "consensus_audit_summary.json").read_text(
            encoding="utf-8"
        )
    )
    selected = points[points["consensus_candidate"]]
    cutoff = summary["view_neighbour_similarity_cutoff"]

    assert (selected["same_seed_structure_similarity"] >= cutoff).all()
    assert (selected["same_seed_electronic_similarity"] >= cutoff).all()
    assert (selected["dual_score"] >= summary["dual_score_threshold"]).all()
    assert (
        selected["same_seed_and_similarity"]
        == selected[
            [
                "same_seed_structure_similarity",
                "same_seed_electronic_similarity",
            ]
        ].min(axis=1)
    ).all()
