import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from analyze_good_materials import joint_rule_performance, numeric_comparison  # noqa: E402


def toy_frame() -> pd.DataFrame:
    n = 100
    zt = np.linspace(0.1, 2.0, n)
    return pd.DataFrame(
        {
            "sample_id": np.arange(n),
            "zt_peak": zt,
            "temperature_peak_K": np.linspace(300, 900, n),
            "abs_seebeck_uV_K": np.linspace(80, 300, n),
            "sigma": np.linspace(1e3, 1e5, n),
            "power_factor_used_mW_mK2": np.linspace(0.2, 6.0, n),
            "kappa_lattice": np.linspace(3.0, 0.1, n),
            "kappa_total": np.linspace(5.0, 0.3, n),
            "material_family": np.where(np.arange(n) % 2, "A", "B"),
        }
    )


def test_numeric_comparison_recovers_expected_directions():
    frame = toy_frame()
    threshold = frame.zt_peak.quantile(0.9)
    result = numeric_comparison(frame, threshold).set_index("feature")
    assert result.loc["power_factor_used_mW_mK2", "auc_probability_good_greater"] > 0.9
    assert result.loc["kappa_total", "auc_probability_good_greater"] < 0.1
    assert result.loc["power_factor_used_mW_mK2", "temperature_strata_used"] >= 1


def test_one_sided_rules_retain_about_ninety_percent_of_good_pairs():
    frame = toy_frame()
    threshold = frame.zt_peak.quantile(0.9)
    numeric = numeric_comparison(frame, threshold)
    rules = joint_rule_performance(frame, threshold, numeric).set_index("rule")
    assert np.isclose(rules.loc["PF soft floor", "high_zt_retention"], 0.9)
    assert np.isclose(rules.loc["kappa_total soft ceiling", "high_zt_retention"], 0.9)
