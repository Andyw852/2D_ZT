"""Phase M: Structure 相似图构建 + k 选择 + fusion sensitivity + 冻结。"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import hellinger_distance, soap_distance, kNN_affinity, graph_qa

root = Path(__file__).resolve().parents[1]
# 载入 SOAP 与 composition
soap_df = pd.read_parquet(root / "features" / "structure" / "geometry_soap_v1.parquet").sort_values("jid").reset_index(drop=True)
comp_df = pd.read_parquet(root / "features" / "structure" / "composition_fraction.parquet").sort_values("jid").reset_index(drop=True)
jids = soap_df["jid"].tolist()
assert jids == comp_df["jid"].tolist()

F = np.array([json.loads(x) for x in comp_df["fraction"]])
d_comp = hellinger_distance(F)
print("d_comp (Hellinger) computed:", d_comp.shape)

# SOAP distances for r_cut 4/6/8
d_geo = {}
for rc in [4, 6, 8]:
    cols = [c for c in soap_df.columns if c.startswith(f"soap{rc}_mean_")]
    X = soap_df[cols].values
    d_geo[rc] = soap_distance(X)
    print(f"d_geo r_cut={rc} computed, range=[{d_geo[rc].min():.4f}, {d_geo[rc].max():.4f}]")

# 归一化到 [0,1]
d_comp_n = d_comp / d_comp.max()
d_geo_n = {rc: d_geo[rc] / d_geo[rc].max() for rc in d_geo}

def knn_overlap(Da, Db, k):
    n = Da.shape[0]
    ka = np.argsort(Da, axis=1)[:, 1:k+1]
    kb = np.argsort(Db, axis=1)[:, 1:k+1]
    return np.mean([len(set(ka[i]) & set(kb[i])) / k for i in range(n)])

# ---- r_cut 稳定性 ----
print("\n=== r_cut 稳定性 (d_geo, k=20) ===")
for a, b in [(4,6),(6,8),(4,8)]:
    print(f"  r_cut {a} vs {b}: kNN(20) overlap = {knn_overlap(d_geo[a], d_geo[b], 20):.4f}")

# ---- baseline: r_cut=6, weight 0.5/0.5, k 扫描 ----
print("\n=== k 扫描 (baseline r_cut=6, w=0.5) ===")
baseline_d = 0.5 * d_geo_n[6] + 0.5 * d_comp_n
for k in [10, 15, 20, 30, 40, 50]:
    W = kNN_affinity(baseline_d, k)
    qa = graph_qa(W)
    print(f"  k={k}: giant={qa['giant_component_fraction']:.4f} comp={qa['n_components']} isolated={qa['isolated_nodes']} mean_deg={qa['mean_degree']}")

# ---- fusion weight 稳定性 (r_cut=6, k=20) ----
print("\n=== fusion weight 稳定性 (r_cut=6, k=20) ===")
d_weights = {}
for w in [0.25, 0.5, 0.75]:
    d_weights[w] = w * d_geo_n[6] + (1 - w) * d_comp_n
for a, b in [(0.25,0.5),(0.5,0.75),(0.25,0.75)]:
    print(f"  w_geo {a} vs {b}: kNN(20) overlap = {knn_overlap(d_weights[a], d_weights[b], 20):.4f}")

# ---- geometry vs composition 邻域差异 (O22 preview) ----
print("\n=== geometry-only vs composition-only (k=20) ===")
print(f"  kNN(20) overlap = {knn_overlap(d_geo[6], d_comp, 20):.4f}")

# 保存 baseline 距离矩阵（后续冻结用）
np.save(root / "data" / "processed" / "d_struct_baseline.npy", baseline_d)
print("\nwrote d_struct_baseline.npy (r_cut=6, w=0.5)")
