import numpy as np
import pytest

from mp_kappaL.te_reference_intersection import (
    empirical_similarity_percentile,
    mean_k_reference_distance,
)


def test_reference_distance_uses_distinct_centroids():
    values = np.array([[0.0, 0.0], [2.0, 0.0]])
    references = np.array([[0.0, 0.0], [1.0, 0.0], [4.0, 0.0]])
    got = mean_k_reference_distance(values, references, k=2)
    np.testing.assert_allclose(got, [0.5, 1.5])


def test_similarity_percentile_is_smaller_distance_better():
    distance = np.array([1.0, 2.0, 3.0, 0.5])
    pool = np.array([True, True, True, False])
    got = empirical_similarity_percentile(distance, pool)
    assert got[3] == 1.0
    assert got[0] > got[1] > got[2]


def test_reference_distance_rejects_invalid_k():
    with pytest.raises(ValueError):
        mean_k_reference_distance(np.zeros((2, 3)), np.zeros((2, 3)), k=3)
