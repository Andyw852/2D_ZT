"""Phase M/N: 构建并冻结所有单视图相似图。"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timezone
from sklearn.preprocessing import RobustScaler
from scipy.spatial.distance import cdist
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import kNN_affinity, graph_qa

root = Path(__file__).resolve().parents[1]
graphs_dir = root / "graphs"
graphs_dir.mkdir(exist_ok=True)
now = datetime.now(timezone.utc).isoformat()

def freeze(name, D, jids, k, meta):
    W = kNN_affinity(D, k)
    qa = graph_qa(W)
    np.savez_compressed(graphs_dir / f"G_{name}.npz", W=W)
    nd = pd.DataFrame({"jid": jids, "degree": np.asarray((W > 0).sum(axis=1)).ravel()})
    nd.to_csv(graphs_dir / f"G_{name}_nodes.csv", index=False)
    full = dict(graph_name=f"G_{name}", k=k, **qa, creation_date=now, random_seed=42,
                knn_tiebreak="fixed-seed(0) deterministic perturbation (<1e-9 * min positive distance): 并列近邻可复现且与 jid 行序无关")
    full.update(meta)
    (graphs_dir / f"G_{name}_metadata.json").write_text(json.dumps(full, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  {name}: N={len(jids)} k={k} giant={qa['giant_component_fraction']:.4f} comp={qa['n_components']} isolated={qa['isolated_nodes']}")
    return W

def scaled_dist(df, cols, jid_col="jid"):
    X = RobustScaler().fit_transform(df[cols].values)
    return cdist(X, X), df[jid_col].tolist()

# ---------- Structure (r_cut=6, w=0.5, k=15) ----------
d_struct = np.load(root / "data" / "processed" / "d_struct_baseline.npy")
soap_df = pd.read_parquet(root / "features" / "structure" / "geometry_soap_v1.parquet").sort_values("jid")
print("=== Structure graph ===")
freeze("structure_v1", d_struct, soap_df["jid"].tolist(), 15,
       {"view":"structure","features":"geometry-only SOAP(r_cut=6,n_max=6,l_max=6,sigma=1.0,periodic=True,mean-pool,L2)+elemental-fraction(Hellinger)",
        "feature_transform":"L2-normalize SOAP + sqrt(2-2K); Hellinger; both normalized to [0,1]",
        "distance_metric":"0.5*d_geo_norm + 0.5*d_comp_norm","kernel":"local-scale Gaussian","symmetrization":"union","normalization":"none"})

# ---------- Electronic graphs ----------
edf = pd.read_parquet(root / "features" / "electronic" / "electronic_features_v1.parquet")
print("=== Electronic graphs ===")
# Eg layer (1103)
eg = edf[["jid","Eg_optb88vdw"]].dropna()
egj = eg["jid"].tolist()
D_eg = np.abs(eg["Eg_optb88vdw"].values[:,None] - eg["Eg_optb88vdw"].values[None,:])
_n_metal = int((eg["Eg_optb88vdw"] == 0).sum())
freeze("Eg_v1", D_eg, egj, 15, {"view":"electronic","features":"Eg_optb88vdw","distance_metric":"|Eg_i-Eg_j|","feature_transform":"raw (0=metal, NOT missing)",
        "note": f"Eg=0 金属 {_n_metal} 个形成零距离等价类：它们之间的 15 近邻由固定种子 tie-break 确定（可复现但金属簇内近邻无物理意义；建议后续用 DOS/有效质量作二阶 tie-break）"})

# electronic-n/p/joint (678)
for name, cols in [("electronic_n_v1", ["Eg_optb88vdw","m_elec_median"]),
                   ("electronic_p_v1", ["Eg_optb88vdw","m_hole_median"]),
                   ("electronic_joint_sensitivity", ["Eg_optb88vdw","m_elec_median","m_hole_median"])]:
    sub = edf.dropna(subset=cols)
    D, jids = scaled_dist(sub, cols)
    freeze(name, D, jids, 15, {"view":"electronic","features":cols,"distance_metric":"Euclidean","feature_transform":"RobustScaler"})

# ---------- Transport graphs ----------
print("=== Transport graphs ===")
V1 = ["S_median","S_MAD","S_sign_fraction","log_sigma_dom_geo","D_sigma","A_sigma_dom"]
for name, f in [("n_transport_v1","n_transport_features_v1.parquet"),
                ("p_transport_v1","p_transport_features_v1.parquet")]:
    tdf = pd.read_parquet(root / "features" / "transport" / f)
    D, jids = scaled_dist(tdf, V1)
    freeze(name, D, jids, 15, {"view":"transport","features":V1,"distance_metric":"Euclidean","feature_transform":"RobustScaler"})

# kappa sensitivity (T2)
T2 = ["S_median","S_MAD","S_range","log_kappa_dom_geo","D_kappa","A_kappa_dom"]
for name, f in [("n_transport_kappa_sensitivity","n_transport_features_candidates.parquet"),
                ("p_transport_kappa_sensitivity","p_transport_features_candidates.parquet")]:
    tdf = pd.read_parquet(root / "features" / "transport" / f)
    D, jids = scaled_dist(tdf, T2)
    freeze(name, D, jids, 15, {"view":"transport_kappa_sensitivity","features":T2,"distance_metric":"Euclidean","feature_transform":"RobustScaler"})

print("\nAll graphs frozen.")
