"""跨视图度量：kNN 并列退化诊断 + 近邻重叠 + 成对距离相关（Step 4 / Step 9）。

Step 4：kNN 带 tiebreak 破平局，并列诊断报告零距离占比与 k/k+1 近邻相等比例。
Step 9：报置换 z 分数 + 效应量 bootstrap 95% CI，不报顶到置换下限的 p 值；
        重叠同时报绝对值与富集倍数；附 Benjamini-Hochberg 多重比较校正。
"""
import hashlib

import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform


def _id_rank(ids):
    """每个 material_id 的稳定 64-bit 散列，用作并列破平局的次排序键。

    旧实现把 UTF-8 字节按小端整数再对 ``2**32`` 取模；对 ``mp-*`` 这类共享
    前四个字节的 ID 会产生大面积碰撞，实际上没有完成稳定破平局。
    """
    return np.array([
        int.from_bytes(
            hashlib.blake2b(str(mid).encode("utf-8"), digest_size=8).digest(),
            "little",
            signed=False,
        )
        for mid in ids
    ], dtype=np.uint64)


def knn_neighbor_matrix(D, k, tiebreak_seed=0, ids=None):
    """返回近邻索引矩阵 (n, k)，第 i 行是 i 的 k 个近邻（不含自身）。

    并列距离的破平局（Step 4）：
    - 提供 ids 时，用 lexsort((id_rank, D[i]))：主键距离、次键按 material 身份，
      可复现且不依赖输入行序（避免 O(n²) 扰动矩阵的内存/时间开销）。
    - 未提供 ids 时，用固定种子扰动破平局。
    """
    D = np.asarray(D, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D 必须是方阵")
    n = D.shape[0]
    if not 0 < k < n:
        raise ValueError(f"k 必须满足 0 < k < n；收到 k={k}, n={n}")
    if ids is not None and len(ids) != n:
        raise ValueError("ids 长度必须与 D 的行数一致")
    nn = np.empty((n, k), dtype=int)
    if ids is not None:
        r = _id_rank(ids)
        for i in range(n):
            order = np.lexsort((r, D[i]))[: k + 1]
            order = order[order != i][: k]
            nn[i] = order
        return nn
    if tiebreak_seed is not None:
        pos = D[D > 0]
        scale = pos.min() if pos.size else 1.0
        D = D + np.random.default_rng(tiebreak_seed).random((n, n)) * scale * 1e-9
        np.fill_diagonal(D, 0.0)
    for i in range(n):
        order = np.argsort(D[i], kind="stable")[: k + 1]
        order = order[order != i][: k]
        nn[i] = order
    return nn


def degeneracy_diagnostics(D, k):
    """Step 4 并列诊断：零非对角元占比 + 每点第 k 与第 k+1 近邻距离相等比例。"""
    n = D.shape[0]
    D = np.asarray(D, dtype=float)
    off = D[~np.eye(n, dtype=bool)]
    zero_frac = float((off == 0.0).mean()) if off.size else 0.0
    tie_frac = 0.0
    if n > k + 1:
        kth = np.partition(D, k, axis=1)[:, k]
        k1th = np.partition(D, k + 1, axis=1)[:, k + 1]
        tie_frac = float((kth == k1th).mean())
    return {"zero_offdiag_frac": zero_frac, "k_k1_tie_frac": tie_frac}


def _row_intersect_count(A, B):
    """逐行 |A[i] ∩ B[i]|，A、B 为 (n, k) 整数近邻索引矩阵。

    -1 是「非共同材料」哨兵：只统计 A 中 >=0 的匹配，避免哨兵与哨兵误匹配。
    """
    n = A.shape[0]
    cnt = np.zeros(n, dtype=int)
    for c in range(A.shape[1]):
        hit = (A[:, c:c + 1] == B) & (A[:, c:c + 1] >= 0)
        cnt += hit.any(axis=1)
    return cnt


def crossview_overlap(NNA, NNB, k, rng, n_perm=1000, n_boot=1000):
    """近邻重叠：绝对重叠 + 置换 z + 富集倍数 + bootstrap 95% CI（Step 9）。

    NNA, NNB 已按共同 material 对齐（第 i 行同一材料）。返回 dict。
    """
    n = NNA.shape[0]
    ov = np.array([np.intersect1d(NNA[i], NNB[i]).size for i in range(n)], dtype=float) / k
    overlap = float(ov.mean())

    # 解析 null mean：E|NNA[i] ∩ NNB[perm[i]]| = Σ_a indegA[a]·indegB[a] / (k·n²)
    # 只统计共同材料（>=0），忽略 -1 哨兵
    indegA = np.bincount(NNA[NNA >= 0].ravel(), minlength=n).astype(float)
    indegB = np.bincount(NNB[NNB >= 0].ravel(), minlength=n).astype(float)
    null_mean = float(indegA @ indegB) / (k * n * n)

    # 置换采样估计 null std（向量化）
    samples = np.empty(n_perm)
    for t in range(n_perm):
        perm = rng.permutation(n)
        samples[t] = _row_intersect_count(NNA, NNB[perm]).mean() / k
    null_std = float(samples.std())
    z = (overlap - null_mean) / (null_std + 1e-12)
    enrich = overlap / (null_mean + 1e-12)

    # bootstrap（按 material_id 重采样）
    boots = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        boots[b] = ov[idx].mean()
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

    return {
        "n": n, "k": k,
        "overlap": overlap,
        "null_mean": null_mean,
        "null_std": null_std,
        "z": z,
        "enrichment": enrich,
        "ci_lo": ci_lo, "ci_hi": ci_hi,
    }


def distance_spearman(DA, DB, rng, n_samp=80_000, n_mantel=50, n_boot=200):
    """成对距离相关：抽样 Spearman + Mantel 置换 z + bootstrap 95% CI（Step 9）。

    DA, DB 为 (n, n) 距离矩阵（已按共同材料对齐）。返回 dict。
    """
    n = DA.shape[0]
    n_pairs = n * (n - 1) // 2
    nsamp = min(n_samp, n_pairs)
    if n_pairs <= n_samp:
        ip, jp = np.triu_indices(n, k=1)
    else:
        # 不构造 O(n²) 的 triu 索引（n=12k 时仅索引就超过 1 GB）。
        ip = rng.randint(0, n, nsamp)
        jp = rng.randint(0, n, nsamp)
        same = ip == jp
        while same.any():
            jp[same] = rng.randint(0, n, int(same.sum()))
            same = ip == jp
    da = DA[ip, jp]
    db = DB[ip, jp]
    rho, _ = stats.spearmanr(da, db)

    # Mantel 置换：行/列同时置换，估计 null mean/std -> z
    nulls = np.empty(n_mantel)
    for t in range(n_mantel):
        perm = rng.permutation(n)
        dbp = DB[perm[ip], perm[jp]]
        nulls[t] = stats.spearmanr(da, dbp).statistic
    null_mean = float(nulls.mean())
    null_std = float(nulls.std())
    z = (rho - null_mean) / (null_std + 1e-12)

    # material-level bootstrap：重采样节点，再从这些节点形成的距离对中抽样。
    # 旧实现直接重采样 pair，把共享端点的 O(n^2) 个距离对误当独立观测，CI 过窄。
    boots = np.empty(n_boot)
    for b in range(n_boot):
        nodes = rng.randint(0, n, n)
        n_pair_boot = min(n_samp, n * (n - 1) // 2)
        bi = rng.randint(0, n, n_pair_boot)
        bj = rng.randint(0, n, n_pair_boot)
        same = bi == bj
        while same.any():
            bj[same] = rng.randint(0, n, int(same.sum()))
            same = bi == bj
        ba = DA[nodes[bi], nodes[bj]]
        bb = DB[nodes[bi], nodes[bj]]
        boots[b] = stats.spearmanr(ba, bb).statistic
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

    return {
        "n": n, "spearman": float(rho),
        "mantel_null_mean": null_mean, "mantel_null_std": null_std,
        "z": float(z), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
    }


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg 校正（Step 9）。输入 p 值列表，返回校正后 p 值数组。"""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / (np.arange(n) + 1.0)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n)
    out[order] = adjusted
    return out


def pairwise_distance_matrix(X, log=False, scale=True):
    """把 (n, d) 特征表转成 (n, n) 欧氏距离矩阵。log 取 log10，scale 用 RobustScaler。"""
    if log:
        X = np.log10(np.clip(X, 1e-12, None))
    if scale and X.shape[1] > 1:
        from sklearn.preprocessing import RobustScaler
        X = RobustScaler().fit_transform(X)
    return squareform(pdist(X)).astype(np.float32)
