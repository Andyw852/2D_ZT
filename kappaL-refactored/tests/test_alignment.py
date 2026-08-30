"""Step 3：验证 load_aligned 能检测并拒绝静默行错位。

覆盖 bug B / bug C 的防御：
- 特征矩阵行数 != row_index 长度 -> 抛异常
- meta 的 material_id 重复 -> 抛异常（bug C）
- 正确时按 material_id 重排特征到 meta 顺序
- inner join 只保留交集，绝不静默错位
"""
import numpy as np
import pandas as pd
import pytest

from mp_kappaL.data_utils import load_aligned
from mp_kappaL.metrics import _id_rank, knn_neighbor_matrix
from mp_kappaL.crossview_analysis import common_pair_neighbors


def _write(tmp_path, ids, X):
    rid = tmp_path / "row_index.npy"
    np.save(rid, np.array(ids))
    feat = tmp_path / "feat.npy"
    np.save(feat, X)
    return str(feat), str(rid)


def test_detects_row_length_mismatch(tmp_path):
    # bug B 型：特征 3 行但 row_index 只有 2 个 id
    feat, rid = _write(tmp_path, ["a", "b"], np.random.RandomState(0).rand(3, 4))
    meta = pd.DataFrame({"material_id": ["a", "b"]})
    with pytest.raises(AssertionError):
        load_aligned({"feat": feat}, meta, row_index_path=rid)


def test_detects_duplicate_meta_id(tmp_path):
    # bug C 型：meta 的 material_id 重复
    feat, rid = _write(tmp_path, ["a", "b", "c"], np.random.RandomState(1).rand(3, 4))
    meta = pd.DataFrame({"material_id": ["a", "a", "c"]})
    with pytest.raises(AssertionError):
        load_aligned({"feat": feat}, meta, row_index_path=rid)


def test_reorders_features_to_meta_order(tmp_path):
    feat, rid = _write(tmp_path, ["a", "b", "c"], np.array([[1.0], [2.0], [3.0]]))
    meta = pd.DataFrame({"material_id": ["c", "a", "b"], "y": [3.0, 1.0, 2.0]})
    m, mats = load_aligned({"feat": feat}, meta, row_index_path=rid)
    assert list(m["material_id"]) == ["c", "a", "b"]
    np.testing.assert_allclose(mats["feat"].ravel(), [3.0, 1.0, 2.0])


def test_inner_join_keeps_intersection_only(tmp_path):
    feat, rid = _write(tmp_path, ["a", "b", "c"], np.array([[1.0], [2.0], [3.0]]))
    # 'd' 不在特征里，'b' 不在 meta 里 -> 交集 {a, c}
    meta = pd.DataFrame({"material_id": ["a", "c", "d"], "y": [1.0, 3.0, 9.0]})
    m, mats = load_aligned({"feat": feat}, meta, row_index_path=rid)
    assert list(m["material_id"]) == ["a", "c"]
    np.testing.assert_allclose(mats["feat"].ravel(), [1.0, 3.0])


def test_material_id_tiebreak_hash_is_not_prefix_colliding():
    ids = [f"mp-{i:08d}" for i in range(100)]
    ranks = _id_rank(ids)
    assert len(np.unique(ranks)) == len(ids)


def test_tied_knn_is_invariant_to_row_order():
    ids = np.array([f"mp-{i}" for i in range(8)])
    # 所有非对角距离都并列；结果只能由 material_id 决定，不能由行号决定。
    D = np.ones((len(ids), len(ids)))
    np.fill_diagonal(D, 0.0)
    base = knn_neighbor_matrix(D, 3, ids=ids)
    base_sets = {ids[i]: set(ids[base[i]]) for i in range(len(ids))}

    perm = np.array([4, 1, 7, 0, 6, 2, 5, 3])
    Dp = D[np.ix_(perm, perm)]
    idp = ids[perm]
    got = knn_neighbor_matrix(Dp, 3, ids=idp)
    got_sets = {idp[i]: set(idp[got[i]]) for i in range(len(idp))}
    assert got_sets == base_sets


def test_crossview_knn_is_recomputed_on_common_cohort():
    ids_a = np.array(["a", "b", "c", "x", "y"])
    xa = np.array([0.0, 10.0, 20.0, 0.1, 9.9])
    da = np.abs(xa[:, None] - xa[None, :])
    ids_b = np.array(["a", "b", "c"])
    xb = np.array([0.0, 10.0, 20.0])
    db = np.abs(xb[:, None] - xb[None, :])
    views = {"A": {"D": da, "ids": ids_a}, "B": {"D": db, "ids": ids_b}}

    common, nna, nnb, _, _ = common_pair_neighbors(views, "A", "B", k=1)
    assert list(common) == ["a", "b", "c"]
    # 两个视图在共同 cohort 上距离完全相同，近邻也必须完全相同。
    np.testing.assert_array_equal(nna, nnb)


if __name__ == "__main__":
    # 无 pytest 时的独立运行（需要 tmp_path，直接用临时目录）
    import tempfile, os
    fns = [test_detects_row_length_mismatch, test_detects_duplicate_meta_id,
           test_reorders_features_to_meta_order, test_inner_join_keeps_intersection_only]
    with tempfile.TemporaryDirectory() as td:
        class Tmp:
            def __init__(self, d): self._d = d
            def __truediv__(self, p): return os.path.join(self._d, p)
        tmp = Tmp(td)
        for fn in fns:
            fn(tmp)
            print(f"PASS {fn.__name__}")
    print(f"all {len(fns)} tests passed")
