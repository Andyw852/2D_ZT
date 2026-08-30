"""共享图工具：kNN + local-scale affinity + QA。"""
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.csgraph import connected_components

def hellinger_distance(F):
    """F: (n, m) 元素分数矩阵。返回 (n, n) Hellinger 距离。"""
    R = np.sqrt(F)
    n = F.shape[0]
    # 用分块避免 O(n^2 m) 过大内存
    D = np.zeros((n, n))
    for i in range(n):
        diff = R[i] - R[i + 1:]
        D[i, i + 1:] = np.sqrt(0.5 * (diff ** 2).sum(axis=1))
    D = D + D.T
    return D

def soap_distance(X):
    """X: (n, d) SOAP mean vectors。L2 normalize 后 d = sqrt(max(0, 2-2K))。"""
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    K = Xn @ Xn.T
    K = np.clip(K, -1, 1)
    return np.sqrt(np.clip(2 - 2 * K, 0, None))

def kNN_affinity(D, k, kernel="local_scale", eps=1e-12, tiebreak_seed=0):
    """D: (n,n) 距离。返回对称化(union)后的 affinity 矩阵 W。

    tiebreak_seed: 并列距离的确定性破平局。距离矩阵里存在大量精确并列
    （如 Eg=0 金属之间 |Eg_i-Eg_j|=0），np.argsort 的稳定排序会让并列近邻
    依赖输入数组顺序（= parquet 的 jid 存储顺序），图结构随 JID 行序漂移。
    这里加一个量级远小于最小正距离的固定种子扰动，使并列近邻选择可复现、
    与输入顺序无关。金属簇内的具体 15 近邻仍无物理意义（见 23 脚本 metadata）。
    """
    n = D.shape[0]
    D = np.asarray(D, dtype=float)
    if tiebreak_seed is not None:
        pos = D[D > 0]
        scale = pos.min() if pos.size else 1.0
        D = D + np.random.default_rng(tiebreak_seed).random((n, n)) * scale * 1e-9
    W = lil_matrix((n, n))
    for i in range(n):
        nn = np.argsort(D[i])[:k + 1]  # 含自身
        nn = nn[nn != i][:k]
        sigma_i = D[i, nn[-1]] if len(nn) else 1.0
        for j in nn:
            sigma_j = D[j, np.argsort(D[j])[:k + 1]][-1] if k > 0 else 1.0
            # 更简单：sigma_j 用 j 的第 k 近邻
            djj = np.sort(D[j])[min(k, n - 1)]
            w = np.exp(-D[i, j] ** 2 / (sigma_i * djj + eps))
            W[i, j] = w
            W[j, i] = w  # union symmetrization
    W = W.tocsr()
    W.setdiag(0)
    return W

def graph_qa(W):
    """W: 对称邻接/affinity。返回连通性统计。"""
    n = W.shape[0]
    n_comp, labels = connected_components(csgraph=W, directed=False)
    comp_sizes = np.bincount(labels)
    giant = comp_sizes.max() / n if n else 0
    degrees = np.asarray((W > 0).sum(axis=1)).ravel()  # unweighted degree
    isolated = int((degrees == 0).sum())
    return {
        "N_nodes": n,
        "N_edges": int((W > 0).nnz) // 2,
        "n_components": n_comp,
        "giant_component_fraction": round(float(giant), 6),
        "isolated_nodes": isolated,
        "mean_degree": round(float(degrees.mean()), 3),
        "median_degree": round(float(np.median(degrees)), 3),
        "min_degree": int(degrees.min()) if n else 0,
        "max_degree": int(degrees.max()) if n else 0,
    }

def neighbor_overlap(Wa, Wb, k, common_idx=None):
    """基于距离/affinity 的 kNN overlap（用 affinity 行排序）。"""
    # 这里输入的是 affinity，直接行排序取 top-k
    n = Wa.shape[0]
    def topk(W, i):
        row = np.asarray(W[i].todense()).ravel() if hasattr(W[i], 'todense') else W[i]
        return set(np.argsort(-row)[:k + 1]) - {i}
    return None  # 由调用方实现
