"""多视图对齐共享工具。"""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix, lil_matrix, diags
from scipy.sparse.csgraph import connected_components
import scipy.sparse.linalg as sla
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import squareform, pdist

def kNN_affinity_from_dist(D, k=15, eps=1e-12):
    n = D.shape[0]
    W = lil_matrix((n, n))
    for i in range(n):
        nn = np.argsort(D[i])[:k + 1]
        nn = nn[nn != i][:k]
        if len(nn) == 0:
            continue
        si = D[i, nn[-1]]
        for j in nn:
            sj = np.sort(D[j])[min(k, n - 1)]
            w = np.exp(-D[i, j] ** 2 / (si * sj + eps))
            W[i, j] = w
            W[j, i] = w
    W = W.tocsr(); W.setdiag(0)
    return W

def load_W(name, gdir):
    d = np.load(gdir / ('G_' + name + '.npz'), allow_pickle=True)
    W = d['W']
    return W.item() if hasattr(W, 'item') else W

def jids_of(name, gdir):
    return pd.read_csv(gdir / ('G_' + name + '_nodes.csv'))['jid'].tolist()

def scale_layer(W):
    """W / mean_strength，使 mean node strength ≈ 1。"""
    deg = np.asarray(W.sum(axis=1)).ravel()
    ms = deg.mean()
    if ms > 1e-12:
        return (W / ms).tocsr(), ms
    return W.tocsr(), ms

def build_supra(layers, consensus_jids, lam):
    """layers: list of (name, sparse W, jids)。
    返回 sparse A, 以及 node 索引信息 dict。"""
    Nc = len(consensus_jids)
    cid = {j: i for i, j in enumerate(consensus_jids)}
    blocks = []  # (global_offset, W_local, jids)
    # 各 view copy 的全局偏移
    offsets = {}
    total = Nc
    for name, W, jids in layers:
        offsets[name] = total
        total += len(jids)
    A = lil_matrix((total, total))
    # block diagonal: view 内部 affinity
    for name, W, jids in layers:
        off = offsets[name]
        A[off:off + W.shape[0], off:off + W.shape[0]] = W
    # identity edges C_i <-> V_i, weight lam / m_i
    for jid in consensus_jids:
        m_i = sum(1 for (name, W, jids) in layers if jid in jids)
        if m_i == 0:
            continue
        ci = cid[jid]
        w = lam / m_i
        for name, W, jids in layers:
            jmap = {j: i for i, j in enumerate(jids)}
            if jid in jmap:
                vi = offsets[name] + jmap[jid]
                A[ci, vi] = w
                A[vi, ci] = w
    return A.tocsr(), offsets, cid

def joint_embedding(A, n_eigs=30):
    deg = np.asarray(A.sum(axis=1)).ravel()
    d_inv_sqrt = 1.0 / np.sqrt(np.maximum(deg, 1e-12))
    D = diags(d_inv_sqrt)
    L = sparse.identity(A.shape[0]) - D @ A @ D
    try:
        ev, V = sla.eigsh(L, k=n_eigs + 1, which='SM', return_eigenvectors=True)
        # 去掉 trivial (最小本征值 ~0 对应常数向量)
        order = np.argsort(ev)
        ev = ev[order]; V = V[:, order]
        return ev, V
    except Exception as e:
        return None, None

def view_coords(V, offsets, name, n_eigs):
    """取出某 view copy 的 joint 坐标（去掉 trivial 后前 n_eigs 个）。"""
    off = offsets[name]
    # V 的第 1 列是 trivial 常数向量，取第 1..n_eigs 列
    return V[off:off + V.shape[0], 1:n_eigs + 1] if False else None
