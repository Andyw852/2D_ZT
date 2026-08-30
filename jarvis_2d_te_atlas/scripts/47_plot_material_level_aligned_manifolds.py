"""Material-level aligned maps for structure, electronic, transport, and kL views.

Every plotted point is a material.  Each view is embedded independently from
its own material-distance matrix and then rigidly Procrustes-aligned to a
rank-distance consensus.  Alignment changes only rotation/reflection/scale;
within-view material geometry is not mixed or distorted.
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


def soap_distance(x: np.ndarray) -> np.ndarray:
    xn = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-15)
    kernel = np.clip(xn @ xn.T, -1.0, 1.0)
    return np.sqrt(np.clip(2.0 - 2.0 * kernel, 0.0, None))


def hellinger_distance(fractions: np.ndarray) -> np.ndarray:
    return squareform(pdist(np.sqrt(fractions), metric="euclidean")) / np.sqrt(2.0)


def table_distance(df: pd.DataFrame, cols: list[str], log_cols: tuple[str, ...] = ()) -> tuple[np.ndarray, list[str]]:
    sub = df.dropna(subset=cols).reset_index(drop=True)
    x = sub[cols].to_numpy(float)
    for col in log_cols:
        j = cols.index(col)
        x[:, j] = np.log10(np.maximum(x[:, j], 1e-15))
    if x.shape[1] > 1:
        x = RobustScaler().fit_transform(x)
    return squareform(pdist(x)), sub["jid"].astype(str).tolist()


def restrict_common(views: dict[str, tuple[np.ndarray, list[str]]]) -> tuple[list[str], dict[str, np.ndarray]]:
    common = sorted(set.intersection(*(set(jids) for _, jids in views.values())))
    matrices = {}
    for name, (distance, jids) in views.items():
        index = {jid: i for i, jid in enumerate(jids)}
        take = [index[jid] for jid in common]
        matrices[name] = distance[np.ix_(take, take)]
    return common, matrices


def rank_normalized_distance(distance: np.ndarray) -> np.ndarray:
    n = distance.shape[0]
    iu = np.triu_indices(n, 1)
    values = rankdata(distance[iu], method="average")
    values = (values - values.min()) / max(values.max() - values.min(), 1e-15)
    out = np.zeros((n, n), dtype=float)
    out[iu] = values
    return out + out.T


def combine_rank_distances(distances: list[np.ndarray]) -> np.ndarray:
    return np.mean([rank_normalized_distance(d) for d in distances], axis=0)


def classical_mds(distance: np.ndarray, dimensions: int = 2) -> np.ndarray:
    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distance ** 2) @ centering
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.clip(values[order][:dimensions], 0.0, None)
    return vectors[:, order][:, :dimensions] * np.sqrt(values)


def normalize_coords(coords: np.ndarray) -> np.ndarray:
    coords = coords - coords.mean(axis=0, keepdims=True)
    return coords / max(np.sqrt(np.mean(np.sum(coords ** 2, axis=1))), 1e-15)


def aligned_material_map(matrices: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, int, float]:
    ranked = {name: rank_normalized_distance(d) for name, d in matrices.items()}
    consensus_distance = np.mean(list(ranked.values()), axis=0)
    consensus_10 = normalize_coords(classical_mds(consensus_distance, dimensions=10))
    consensus_2 = normalize_coords(consensus_10[:, :2])
    aligned = {}
    for name, distance in ranked.items():
        coords = normalize_coords(classical_mds(distance, dimensions=2))
        rotation, _ = orthogonal_procrustes(coords, consensus_2)
        aligned[name] = normalize_coords(coords @ rotation)
    aligned["Consensus"] = consensus_2

    scores = []
    fits = {}
    for k in range(4, 9):
        model = KMeans(n_clusters=k, n_init=50, random_state=SEED).fit(consensus_10)
        score = silhouette_score(consensus_10, model.labels_)
        scores.append((score, k))
        fits[k] = model
    best_score, best_k = max(scores)
    raw = fits[best_k].labels_
    means = pd.DataFrame({"label": raw, "x": consensus_2[:, 0], "y": consensus_2[:, 1]}).groupby("label").mean()
    order = means.sort_values(["x", "y"]).index.tolist()
    remap = {old: new + 1 for new, old in enumerate(order)}
    labels = np.asarray([remap[x] for x in raw], dtype=int)
    return aligned, consensus_10, labels, best_k, float(best_score)


def medoids(consensus_10: np.ndarray, labels: np.ndarray) -> set[int]:
    selected = set()
    for label in sorted(set(labels)):
        take = np.where(labels == label)[0]
        center = consensus_10[take].mean(axis=0)
        selected.add(int(take[np.argmin(np.linalg.norm(consensus_10[take] - center, axis=1))]))
    return selected


def load_2d_atlas() -> tuple[list[str], dict[str, np.ndarray], pd.DataFrame]:
    structure_nodes = pd.read_csv(ROOT / "graphs/G_structure_v1_nodes.csv").sort_values("jid")
    structure_distance = np.load(ROOT / "data/processed/d_struct_baseline.npy")
    electronic = pd.read_parquet(ROOT / "features/electronic/electronic_features_v1.parquet")
    views: dict[str, tuple[np.ndarray, list[str]]] = {
        "Structure": (structure_distance, structure_nodes["jid"].astype(str).tolist()),
        "Electronic n": table_distance(electronic, ["Eg_optb88vdw", "m_elec_median"]),
        "Electronic p": table_distance(electronic, ["Eg_optb88vdw", "m_hole_median"]),
    }
    transport_cols = ["S_median", "S_MAD", "S_sign_fraction", "log_sigma_dom_geo", "D_sigma", "A_sigma_dom"]
    for carrier in ("n", "p"):
        table = pd.read_parquet(ROOT / f"features/transport/{carrier}_transport_features_v1.parquet")
        views[f"Transport {carrier}"] = table_distance(table, transport_cols)
    jids, raw = restrict_common(views)
    matrices = {
        "Structure": raw["Structure"],
        "Electronic": combine_rank_distances([raw["Electronic n"], raw["Electronic p"]]),
        "Transport": combine_rank_distances([raw["Transport n"], raw["Transport p"]]),
    }
    structures = pd.read_parquet(ROOT / "data/processed/standardized_2d_structures.parquet")
    meta = structures[["jid", "formula"]].drop_duplicates("jid").set_index("jid").reindex(jids).reset_index()
    n_pf = pd.read_parquet(ROOT / "features/transport/n_transport_features_v1.parquet")[["jid", "PF_mean"]].rename(columns={"PF_mean": "PF_n"})
    p_pf = pd.read_parquet(ROOT / "features/transport/p_transport_features_v1.parquet")[["jid", "PF_mean"]].rename(columns={"PF_mean": "PF_p"})
    meta = meta.merge(n_pf, on="jid", how="left").merge(p_pf, on="jid", how="left")
    return jids, matrices, meta


def load_kl_subset() -> tuple[list[str], dict[str, np.ndarray], pd.DataFrame]:
    df = pd.read_parquet(ROOT / "features/kl_verify/kl_views.parquet").reset_index(drop=True)
    geometry = soap_distance(np.load(ROOT / "data/processed/kl_soap_geo.npy"))
    composition = hellinger_distance(np.load(ROOT / "data/processed/kl_comp_frac.npy"))
    structure = 0.5 * geometry / max(geometry.max(), 1e-15) + 0.5 * composition / max(composition.max(), 1e-15)
    views = {
        "Structure": (structure, df["jid"].astype(str).tolist()),
        "Electronic": table_distance(df, ["Eg_opt", "m_elec", "m_hole"]),
        "Elastic": table_distance(df, ["B_kv", "G_gv"]),
        "Lattice kL": table_distance(df, ["kL_300"], log_cols=("kL_300",)),
    }
    jids, matrices = restrict_common(views)
    meta = df[["jid", "formula", "kL_300"]].set_index("jid").reindex(jids).reset_index()
    dual = pd.read_csv(ROOT / "data/processed/pf_kL_dual_channel_intersection.csv")
    pf = dual.groupby("jid", as_index=False).agg(PF_max=("PF_jarvis", "max"), dual_candidate=("top20_intersection", "max"))
    meta = meta.merge(pf, on="jid", how="left")
    return jids, matrices, meta


def build_rows(dataset: str, jids: list[str], matrices: dict[str, np.ndarray], meta: pd.DataFrame) -> tuple[list[dict], dict]:
    aligned, consensus_10, labels, best_k, best_score = aligned_material_map(matrices)
    representative = medoids(consensus_10, labels)
    rows = []
    for i, jid in enumerate(jids):
        row = {
            "dataset": dataset,
            "jid": jid,
            "formula": str(meta.iloc[i].get("formula", jid)),
            "cluster": int(labels[i]),
            "representative": i in representative,
        }
        for name, coords in aligned.items():
            key = name.lower().replace(" ", "_")
            row[f"{key}_x"] = round(float(coords[i, 0]), 6)
            row[f"{key}_y"] = round(float(coords[i, 1]), 6)
        for col in ("PF_n", "PF_p", "PF_max", "kL_300", "dual_candidate"):
            if col in meta.columns:
                value = meta.iloc[i][col]
                row[col] = None if pd.isna(value) else (bool(value) if col == "dual_candidate" else round(float(value), 6))
        rows.append(row)
    info = {"dataset": dataset, "n": len(rows), "views": list(aligned), "k": best_k, "silhouette": best_score}
    return rows, info


def draw_static(all_rows: list[dict], infos: list[dict]) -> None:
    palette = plt.get_cmap("tab10")
    for info in infos:
        dataset = info["dataset"]
        frame = pd.DataFrame([r for r in all_rows if r["dataset"] == dataset])
        views = info["views"]
        ncols = len(views)
        fig, axes = plt.subplots(1, ncols, figsize=(4.4 * ncols, 4.5), constrained_layout=True)
        for ax, view in zip(np.atleast_1d(axes), views):
            key = view.lower().replace(" ", "_")
            for cluster in sorted(frame.cluster.unique()):
                sub = frame[frame.cluster == cluster]
                ax.scatter(sub[f"{key}_x"], sub[f"{key}_y"], s=13, alpha=0.62, color=palette(cluster - 1), linewidths=0)
            for _, row in frame[frame.representative].iterrows():
                ax.annotate(row.formula, (row[f"{key}_x"], row[f"{key}_y"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
            ax.set_title(f"{view}: every point is a material", loc="left", fontsize=10)
            ax.set_xlabel("Aligned dimension 1")
            ax.set_ylabel("Aligned dimension 2")
            ax.grid(alpha=0.15, linewidth=0.5)
        fig.suptitle(f"Material-level aligned feature manifolds — {dataset} (N={len(frame)})", fontsize=14)
        out = "material_level_2d_manifolds" if dataset == "2d_atlas" else "material_level_kl_manifolds"
        fig.savefig(ROOT / f"figures/{out}.png", dpi=240)
        fig.savefig(ROOT / f"figures/{out}.pdf")
        plt.close(fig)


def main() -> None:
    all_rows, infos = [], []
    for dataset, loader in (("2d_atlas", load_2d_atlas), ("kl_subset", load_kl_subset)):
        jids, matrices, meta = loader()
        rows, info = build_rows(dataset, jids, matrices, meta)
        all_rows.extend(rows)
        infos.append(info)
    payload = {
        "method": "view-specific rank-distance MDS; rigid Procrustes alignment to equal-view consensus",
        "datasets": infos,
        "materials": all_rows,
    }
    (ROOT / "data/processed/material_level_aligned_manifolds.json").write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    pd.DataFrame(all_rows).to_parquet(ROOT / "manifolds/material_level_aligned_manifolds.parquet", index=False)
    draw_static(all_rows, infos)
    print(pd.DataFrame(infos).to_string(index=False))
    print("Wrote material-level aligned maps, JSON, parquet, PNG, and PDF.")


if __name__ == "__main__":
    main()
