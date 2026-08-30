"""Build n/p-separated, material-level manifolds with performance overlays.

Every point is a material.  Structure, carrier-specific electronic, and
carrier-specific transport distances are embedded independently, then rigidly
aligned to a carrier-specific equal-view consensus.  The matched experimental
kappa_L subset adds a fourth lattice-kappa view.  PF is an external colour/
label variable and is not used to construct the full-data manifold.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import RobustScaler


ROOT = Path(__file__).resolve().parents[1]
SEED = 42
TRANSPORT_COLS = [
    "S_median",
    "S_MAD",
    "S_sign_fraction",
    "log_sigma_dom_geo",
    "D_sigma",
    "A_sigma_dom",
]


def table_distance(df: pd.DataFrame, cols: list[str], log_cols: tuple[str, ...] = ()):
    sub = df.dropna(subset=cols).reset_index(drop=True)
    x = sub[cols].to_numpy(float).copy()
    for col in log_cols:
        j = cols.index(col)
        x[:, j] = np.log10(np.maximum(x[:, j], 1e-15))
    if x.shape[1] > 1:
        x = RobustScaler().fit_transform(x)
    return squareform(pdist(x)), sub["jid"].astype(str).tolist()


def soap_distance(x: np.ndarray) -> np.ndarray:
    x = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-15)
    kernel = np.clip(x @ x.T, -1, 1)
    return np.sqrt(np.clip(2 - 2 * kernel, 0, None))


def hellinger_distance(x: np.ndarray) -> np.ndarray:
    return squareform(pdist(np.sqrt(x), metric="euclidean")) / np.sqrt(2)


def restrict_common(views):
    common = sorted(set.intersection(*(set(jids) for _, jids in views.values())))
    matrices = {}
    for name, (distance, jids) in views.items():
        index = {jid: i for i, jid in enumerate(jids)}
        take = [index[jid] for jid in common]
        matrices[name] = distance[np.ix_(take, take)]
    return common, matrices


def rank_distance(distance: np.ndarray) -> np.ndarray:
    n = len(distance)
    iu = np.triu_indices(n, 1)
    values = rankdata(distance[iu], method="average")
    values = (values - values.min()) / max(values.max() - values.min(), 1e-15)
    out = np.zeros((n, n), float)
    out[iu] = values
    return out + out.T


def classical_mds(distance: np.ndarray, dimensions: int = 2) -> np.ndarray:
    n = len(distance)
    h = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * h @ (distance**2) @ h
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.clip(values[order][:dimensions], 0, None)
    return vectors[:, order][:, :dimensions] * np.sqrt(values)


def normalize(x: np.ndarray) -> np.ndarray:
    x = x - x.mean(axis=0, keepdims=True)
    return x / max(np.sqrt(np.mean(np.sum(x**2, axis=1))), 1e-15)


def embed(matrices: dict[str, np.ndarray]):
    ranked = {name: rank_distance(d) for name, d in matrices.items()}
    consensus_distance = np.mean(list(ranked.values()), axis=0)
    consensus_10 = normalize(classical_mds(consensus_distance, 10))
    consensus_2 = normalize(consensus_10[:, :2])
    aligned = {}
    for name, distance in ranked.items():
        xy = normalize(classical_mds(distance, 2))
        rotation, _ = orthogonal_procrustes(xy, consensus_2)
        aligned[name] = normalize(xy @ rotation)
    aligned["Consensus"] = consensus_2

    model = KMeans(n_clusters=4, n_init=50, random_state=SEED).fit(consensus_10)
    labels0 = model.labels_
    centers = pd.DataFrame({"label": labels0, "x": consensus_2[:, 0], "y": consensus_2[:, 1]}).groupby("label").mean()
    order = centers.sort_values(["x", "y"]).index.tolist()
    remap = {old: new + 1 for new, old in enumerate(order)}
    labels = np.asarray([remap[x] for x in labels0], int)
    return aligned, labels, float(silhouette_score(consensus_10, labels0))


def base_views(carrier: str):
    nodes = pd.read_csv(ROOT / "graphs/G_structure_v1_nodes.csv").sort_values("jid")
    structure = np.load(ROOT / "data/processed/d_struct_baseline.npy")
    electronic = pd.read_parquet(ROOT / "features/electronic/electronic_features_v1.parquet")
    mass = "m_elec_median" if carrier == "n" else "m_hole_median"
    transport = pd.read_parquet(ROOT / f"features/transport/{carrier}_transport_features_v1.parquet")
    views = {
        "Structure": (structure, nodes["jid"].astype(str).tolist()),
        "Electronic": table_distance(electronic, ["Eg_optb88vdw", mass]),
        "Transport": table_distance(transport, TRANSPORT_COLS),
    }
    jids, matrices = restrict_common(views)
    structures = pd.read_parquet(ROOT / "data/processed/standardized_2d_structures.parquet")
    meta = structures[["jid", "formula"]].drop_duplicates("jid")
    meta = meta.merge(transport[["jid", "PF_mean"]], on="jid", how="left").set_index("jid").reindex(jids).reset_index()
    return jids, matrices, meta


def matched_views(carrier: str):
    df = pd.read_parquet(ROOT / "features/kl_verify/kl_views.parquet").reset_index(drop=True)
    geometry = soap_distance(np.load(ROOT / "data/processed/kl_soap_geo.npy"))
    composition = hellinger_distance(np.load(ROOT / "data/processed/kl_comp_frac.npy"))
    structure = 0.5 * geometry / max(geometry.max(), 1e-15) + 0.5 * composition / max(composition.max(), 1e-15)
    mass = "m_elec" if carrier == "n" else "m_hole"
    views = {
        "Structure": (structure, df["jid"].astype(str).tolist()),
        "Electronic": table_distance(df, ["Eg_opt", mass]),
        "Lattice kL": table_distance(df, ["kL_300"], log_cols=("kL_300",)),
    }
    keep, out = restrict_common(views)
    dual = pd.read_csv(ROOT / "data/processed/pf_kL_dual_channel_intersection.csv")
    dual = dual[dual["carrier"] == carrier].drop_duplicates("jid").set_index("jid")
    keep = [jid for jid in keep if jid in dual.index]
    source_jids = restrict_common(views)[0]
    index = {jid: i for i, jid in enumerate(source_jids)}
    take = [index[jid] for jid in keep]
    out = {name: d[np.ix_(take, take)] for name, d in out.items()}
    m = df[["jid", "formula"]].drop_duplicates("jid").set_index("jid").reindex(keep).reset_index()
    extra_cols = ["PF_jarvis", "kL_exp_300K", "PF_percentile", "low_kL_percentile", "top20_intersection", "pareto"]
    extra = dual.loc[keep, extra_cols].reset_index()
    m = m.merge(extra, on="jid", how="left").rename(columns={"PF_jarvis": "PF_mean"})
    return keep, out, m


def build_frame(dataset: str, carrier: str, jids: list[str], matrices: dict[str, np.ndarray], meta: pd.DataFrame):
    aligned, clusters, silhouette = embed(matrices)
    frame = meta.copy()
    frame["dataset"] = dataset
    frame["carrier"] = carrier
    frame["cluster"] = clusters
    for view, xy in aligned.items():
        key = view.lower().replace(" ", "_")
        frame[f"{key}_x"] = xy[:, 0]
        frame[f"{key}_y"] = xy[:, 1]
    frame["PF_percentile"] = frame.get("PF_percentile", frame["PF_mean"].rank(pct=True))
    if "kL_exp_300K" not in frame:
        frame["kL_exp_300K"] = np.nan
    if "low_kL_percentile" not in frame:
        frame["low_kL_percentile"] = np.nan
    if "top20_intersection" not in frame:
        frame["top20_intersection"] = False
    if "pareto" not in frame:
        frame["pareto"] = False
    frame["top20_intersection"] = frame["top20_intersection"].fillna(False).astype(bool)
    frame["pareto"] = frame["pareto"].fillna(False).astype(bool)
    frame["label"] = False
    if dataset == "full":
        for cluster in sorted(frame["cluster"].unique()):
            top = frame[frame["cluster"] == cluster].nlargest(2, "PF_mean").index
            frame.loc[top, "label"] = True
    else:
        frame["label"] = frame["top20_intersection"] | frame["pareto"]
    info = {"dataset": dataset, "carrier": carrier, "n": len(frame), "views": list(aligned), "silhouette": silhouette}
    return frame, info


def align_p_to_n(frames: dict[str, pd.DataFrame], dataset: str):
    n = frames[f"{dataset}_n"]
    p = frames[f"{dataset}_p"]
    common = sorted(set(n["jid"]) & set(p["jid"]))
    ni = n.set_index("jid").loc[common]
    pi = p.set_index("jid").loc[common]
    a = pi[["consensus_x", "consensus_y"]].to_numpy(float)
    b = ni[["consensus_x", "consensus_y"]].to_numpy(float)
    rotation, _ = orthogonal_procrustes(a, b)
    for key in [c[:-2] for c in p.columns if c.endswith("_x")]:
        xy = p[[f"{key}_x", f"{key}_y"]].to_numpy(float) @ rotation
        p[[f"{key}_x", f"{key}_y"]] = xy


def plot_static(frames: dict[str, pd.DataFrame]):
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    views = ["Structure", "Electronic", "Transport", "Consensus"]
    for row, carrier in enumerate(("n", "p")):
        frame = frames[f"full_{carrier}"]
        color = np.log10(np.maximum(frame["PF_mean"], 1e-9))
        for ax, view in zip(axes[row], views):
            key = view.lower()
            points = ax.scatter(frame[f"{key}_x"], frame[f"{key}_y"], c=color, s=14, cmap="viridis", alpha=0.72)
            for _, r in frame[frame["label"]].iterrows():
                ax.annotate(r["formula"], (r[f"{key}_x"], r[f"{key}_y"]), fontsize=7)
            ax.set_title(f"{carrier}-type · {view}", loc="left")
            ax.set_xlabel("aligned dimension 1")
            ax.set_ylabel("aligned dimension 2")
    fig.colorbar(points, ax=axes, label="log10(PF_mean)", shrink=0.7)
    fig.savefig(ROOT / "figures/np_separated_performance_manifolds.png", dpi=220)
    fig.savefig(ROOT / "figures/np_separated_performance_manifolds.pdf")
    plt.close(fig)


def main():
    frames = {}
    infos = []
    for carrier in ("n", "p"):
        jids, matrices, meta = base_views(carrier)
        full, info = build_frame("full", carrier, jids, matrices, meta)
        frames[f"full_{carrier}"] = full
        infos.append(info)
        kj, km, kmeta = matched_views(carrier)
        matched, info = build_frame("kl", carrier, kj, km, kmeta)
        frames[f"kl_{carrier}"] = matched
        infos.append(info)
    align_p_to_n(frames, "full")
    align_p_to_n(frames, "kl")
    all_rows = pd.concat(frames.values(), ignore_index=True)
    all_rows.to_parquet(ROOT / "manifolds/np_separated_performance_manifolds.parquet", index=False)
    payload = {
        "method": "carrier-specific rank-distance MDS; rigid Procrustes alignment; PF is an external label",
        "datasets": infos,
        "materials": all_rows.where(pd.notna(all_rows), None).to_dict(orient="records"),
    }
    (ROOT / "data/processed/np_separated_performance_manifolds.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_static(frames)
    for info in infos:
        print(info)


if __name__ == "__main__":
    main()
