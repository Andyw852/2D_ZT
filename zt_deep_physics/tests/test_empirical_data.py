from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "empirical"))

from build_empirical_atlas import (  # noqa: E402
    clean_curve,
    parse_grain_size_um,
    parse_relative_density,
)


def test_relative_density_parser_handles_threshold_and_range():
    assert parse_relative_density(">=95%", "") == 95
    assert parse_relative_density("%", "between 85% and 90%") == 87.5


def test_grain_size_parser_converts_units_and_ranges():
    assert np.isclose(parse_grain_size_um("100nm", "100 to 500 nm"), 0.3)
    assert np.isclose(parse_grain_size_um("1um", "0.5 to 1.5 micrometer"), 1.0)


def test_curve_cleaning_removes_bad_temperature_and_property_values():
    x = np.array([50.0, 300.0, 600.0, 1700.0])
    y = np.array([0.5, 0.8, 7.0, 1.0])
    clean_x, clean_y = clean_curve("zt", x, y)
    assert np.array_equal(clean_x, np.array([300.0]))
    assert np.array_equal(clean_y, np.array([0.8]))

