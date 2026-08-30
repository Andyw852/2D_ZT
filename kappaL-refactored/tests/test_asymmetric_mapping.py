import numpy as np

from mp_kappaL.asymmetric_mapping import composition_features
from mp_kappaL.unified_chemical_space import _unit_inertia


def test_composition_features_are_stoichiometry_scale_invariant():
    x_small, names, group = composition_features(["Al", "Al", "O", "O", "O"])
    x_large, names_large, group_large = composition_features(
        ["Al", "Al", "O", "O", "O"] * 4
    )
    assert x_small.shape == (162,)
    assert len(names) == 162
    assert names == names_large
    assert group == group_large == "Al-O"
    np.testing.assert_allclose(x_small, x_large, rtol=0, atol=1e-7)


def test_composition_features_keep_atomic_fractions():
    x, _, _ = composition_features(["Al", "Al", "O", "O", "O"])
    assert np.isclose(x[12], 0.4)  # Al, Z=13
    assert np.isclose(x[7], 0.6)   # O, Z=8
    assert np.isclose(x[:118].sum(), 1.0)


def test_superblocks_are_normalized_to_equal_total_inertia():
    a = np.arange(30, dtype=float).reshape(10, 3)
    b = np.arange(50, dtype=float).reshape(10, 5) ** 2
    a_norm = _unit_inertia(a)
    b_norm = _unit_inertia(b)
    assert np.allclose(a_norm.mean(axis=0), 0)
    assert np.allclose(b_norm.mean(axis=0), 0)
    assert np.isclose(np.linalg.norm(a_norm, ord="fro"), 1)
    assert np.isclose(np.linalg.norm(b_norm, ord="fro"), 1)
