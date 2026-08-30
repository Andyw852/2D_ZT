"""Step 0 冒烟测试：验证共享图工具 hellinger_distance / soap_distance 的基本性质。

性质：对称性、非负性、对角为零、手算小例子。
"""
import numpy as np
from graph_utils import hellinger_distance, soap_distance


def test_hellinger_small_hand_computed():
    # Hellinger 距离 = sqrt(0.5 * sum((sqrt p - sqrt q)^2))
    F = np.array([
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
    ])
    D = hellinger_distance(F)
    assert D.shape == (2, 2)
    assert np.isclose(D[0, 1], np.sqrt(0.5))
    assert np.isclose(D[1, 0], np.sqrt(0.5))


def test_hellinger_diagonal_zero():
    F = np.array([[0.3, 0.7], [0.3, 0.7]])
    D = hellinger_distance(F)
    assert np.allclose(D[0, 1], 0.0)
    assert np.allclose(np.diag(D), 0.0)


def test_hellinger_symmetric_nonnegative():
    rng = np.random.RandomState(0)
    F = rng.dirichlet(np.ones(5), size=20)
    D = hellinger_distance(F)
    assert D.shape == (20, 20)
    assert np.allclose(D, D.T)
    assert (D >= -1e-12).all()


def test_soap_diagonal_zero():
    rng = np.random.RandomState(1)
    X = rng.rand(10, 8)
    D = soap_distance(X)
    assert D.shape == (10, 10)
    assert np.allclose(np.diag(D), 0.0, atol=1e-6)


def test_soap_identical_vectors_zero():
    X = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
    D = soap_distance(X)
    assert np.allclose(D[0, 1], 0.0, atol=1e-6)


def test_soap_symmetric_nonnegative():
    rng = np.random.RandomState(2)
    X = rng.rand(15, 12)
    D = soap_distance(X)
    assert np.allclose(D, D.T)
    assert (D >= -1e-12).all()


if __name__ == "__main__":
    # 无 pytest 时的独立运行入口
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"all {len(fns)} tests passed")
