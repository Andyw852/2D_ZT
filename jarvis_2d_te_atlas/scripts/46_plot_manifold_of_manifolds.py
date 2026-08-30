"""Build a defensible "manifold of manifolds" distance map.

Each point cloud represents one complete feature view, not one material class.
For a common set of materials we vectorize the within-view material distance
matrix, rank its entries, and define inter-view distance as

    d(A, B) = sqrt(2 * (1 - Spearman(vec(D_A), vec(D_B)))).

Classical MDS embeds those inter-view distances.  Repeated 80% material
subsamples show sampling uncertainty around each view location.  The axes have
no standalone physical meaning; only separations and neighborhoods do.
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
from sklearn.preprocessing import RobustScaler


ROOT = Path(__file__).resolve().parents[1]
SEED = 42
N_BOOT = 300
BOOT_FRACTION = 0.80


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


def common_view_matrices(views: dict[str, tuple[np.ndarray, list[str]]]) -> tuple[list[str], dict[str, np.ndarray]]:
    common = sorted(set.intersection(*(set(jids) for _, jids in views.values())))
    out: dict[str, np.ndarray] = {}
    for name, (distance, jids) in views.items():
        index = {jid: i for i, jid in enumerate(jids)}
        take = [index[jid] for jid in common]
        out[name] = distance[np.ix_(take, take)]
    return common, out


def correlation_geometry(matrices: dict[str, np.ndarray], take: np.ndarray | None = None) -> tuple[list[str], np.ndarray, np.ndarray]:
    names = list(matrices)
    if take is None:
        n = next(iter(matrices.values())).shape[0]
        take = np.arange(n)
    iu = np.triu_indices(len(take), 1)
    signatures = []
    for name in names:
        d = matrices[name][np.ix_(take, take)][iu]
        ranked = rankdata(d, method="average")
        ranked = (ranked - ranked.mean()) / max(ranked.std(), 1e-15)
        signatures.append(ranked)
    signatures = np.asarray(signatures)
    corr = np.clip(np.corrcoef(signatures), -1.0, 1.0)
    distance = np.sqrt(np.clip(2.0 * (1.0 - corr), 0.0, None))
    np.fill_diagonal(distance, 0.0)
    return names, corr, distance


def classical_mds(distance: np.ndarray, dimensions: int = 2) -> np.ndarray:
    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distance ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.clip(eigenvalues[order][:dimensions], 0.0, None)
    eigenvectors = eigenvectors[:, order][:, :dimensions]
    return eigenvectors * np.sqrt(eigenvalues)


def orient(coords: np.ndarray, names: list[str], left: str, right: list[str]) -> np.ndarray:
    coords = coords - coords.mean(axis=0, keepdims=True)
    left_point = coords[names.index(left)]
    right_point = np.mean([coords[names.index(x)] for x in right], axis=0)
    vector = right_point - left_point
    angle = np.arctan2(vector[1], vector[0])
    rotation = np.array([[np.cos(-angle), -np.sin(-angle)], [np.sin(-angle), np.cos(-angle)]])
    coords = coords @ rotation.T
    if coords[names.index(left), 0] > np.mean([coords[names.index(x), 0] for x in right]):
        coords[:, 0] *= -1
    return coords


def bootstrap_map(
    matrices: dict[str, np.ndarray],
    left: str,
    right: list[str],
    rng: np.random.Generator,
) -> dict:
    names, corr, distance = correlation_geometry(matrices)
    base = orient(classical_mds(distance), names, left, right)
    n = next(iter(matrices.values())).shape[0]
    sample_n = max(20, int(round(n * BOOT_FRACTION)))
    boot = np.zeros((N_BOOT, len(names), 2), dtype=float)
    for b in range(N_BOOT):
        take = np.sort(rng.choice(n, sample_n, replace=False))
        _, _, boot_distance = correlation_geometry(matrices, take)
        coords = classical_mds(boot_distance)
        coords -= coords.mean(axis=0, keepdims=True)
        rotation, _ = orthogonal_procrustes(coords, base)
        boot[b] = coords @ rotation
    return {
        "names": names,
        "corr": corr,
        "distance": distance,
        "coords": base,
        "boot": boot,
        "n_materials": n,
        "sample_n": sample_n,
    }


def load_main_views() -> tuple[list[str], dict[str, np.ndarray]]:
    structure_nodes = pd.read_csv(ROOT / "graphs/G_structure_v1_nodes.csv").sort_values("jid")
    structure_distance = np.load(ROOT / "data/processed/d_struct_baseline.npy")
    assert len(structure_nodes) == structure_distance.shape[0]
    views: dict[str, tuple[np.ndarray, list[str]]] = {
        "Structure": (structure_distance, structure_nodes["jid"].astype(str).tolist())
    }
    electronic = pd.read_parquet(ROOT / "features/electronic/electronic_features_v1.parquet")
    views["Band gap"] = table_distance(electronic, ["Eg_optb88vdw"])
    views["Electronic n"] = table_distance(electronic, ["Eg_optb88vdw", "m_elec_median"])
    views["Electronic p"] = table_distance(electronic, ["Eg_optb88vdw", "m_hole_median"])
    transport_cols = ["S_median", "S_MAD", "S_sign_fraction", "log_sigma_dom_geo", "D_sigma", "A_sigma_dom"]
    for carrier in ("n", "p"):
        table = pd.read_parquet(ROOT / f"features/transport/{carrier}_transport_features_v1.parquet")
        views[f"Transport {carrier}"] = table_distance(table, transport_cols)
    return common_view_matrices(views)


def load_kl_views() -> tuple[list[str], dict[str, np.ndarray]]:
    df = pd.read_parquet(ROOT / "features/kl_verify/kl_views.parquet").reset_index(drop=True)
    soap = np.load(ROOT / "data/processed/kl_soap_geo.npy")
    comp = np.load(ROOT / "data/processed/kl_comp_frac.npy")
    geometry = soap_distance(soap)
    composition = hellinger_distance(comp)
    structure = 0.5 * geometry / max(geometry.max(), 1e-15) + 0.5 * composition / max(composition.max(), 1e-15)
    views: dict[str, tuple[np.ndarray, list[str]]] = {
        "Structure": (structure, df["jid"].astype(str).tolist()),
        "Band gap": table_distance(df, ["Eg_opt"]),
        "Electronic": table_distance(df, ["Eg_opt", "m_elec", "m_hole"]),
        "Elastic": table_distance(df, ["B_kv", "G_gv"]),
        "Lattice kL": table_distance(df, ["kL_300"], log_cols=("kL_300",)),
    }
    return common_view_matrices(views)


def rows_for_export(panel: str, result: dict) -> tuple[list[dict], list[dict]]:
    points = []
    for i, name in enumerate(result["names"]):
        points.append({
            "panel": panel,
            "view": name,
            "kind": "centroid",
            "x": float(result["coords"][i, 0]),
            "y": float(result["coords"][i, 1]),
            "n_materials": int(result["n_materials"]),
        })
        for b in range(N_BOOT):
            points.append({
                "panel": panel,
                "view": name,
                "kind": "bootstrap",
                "x": float(result["boot"][b, i, 0]),
                "y": float(result["boot"][b, i, 1]),
                "n_materials": int(result["sample_n"]),
            })
    pairs = []
    for i, a in enumerate(result["names"]):
        for j in range(i + 1, len(result["names"])):
            b = result["names"][j]
            pairs.append({
                "panel": panel,
                "view_a": a,
                "view_b": b,
                "n_common": int(result["n_materials"]),
                "spearman": float(result["corr"][i, j]),
                "distance": float(result["distance"][i, j]),
            })
    return points, pairs


def draw_static(results: dict[str, dict]) -> None:
    palette = {
        "Structure": "#7e57c2",
        "Band gap": "#f39c12",
        "Electronic n": "#16a085",
        "Electronic p": "#2e86c1",
        "Transport n": "#c0392b",
        "Transport p": "#d35400",
        "Electronic": "#168c78",
        "Elastic": "#3f7fbf",
        "Lattice kL": "#c13b56",
    }
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), constrained_layout=True)
    titles = {
        "2d_atlas": "2D atlas: geometry of six feature manifolds",
        "kl_validation": "kL subset: electronic–phonon manifold separation",
    }
    offsets = {
        "2d_atlas": {
            "Structure": (8, 8), "Band gap": (-8, 10), "Electronic n": (8, -12),
            "Electronic p": (-8, 10), "Transport n": (8, -12), "Transport p": (-8, 10),
        },
        "kl_validation": {
            "Structure": (-8, 10), "Band gap": (-8, -14), "Electronic": (8, 10),
            "Elastic": (8, 10), "Lattice kL": (8, 10),
        },
    }
    for label, ax in zip(("2d_atlas", "kl_validation"), axes):
        result = results[label]
        all_boot = result["boot"].reshape(-1, 2)
        ax.scatter(all_boot[:, 0], all_boot[:, 1], s=5, c="0.75", alpha=0.10, linewidths=0)
        for i, name in enumerate(result["names"]):
            cloud = result["boot"][:, i, :]
            color = palette[name]
            ax.scatter(cloud[:, 0], cloud[:, 1], s=8, color=color, alpha=0.16, linewidths=0)
            x, y = result["coords"][i]
            ax.scatter([x], [y], s=80, color=color, edgecolor="white", linewidth=1.0, zorder=4)
            dx, dy = offsets[label][name]
            ax.annotate(
                name,
                (x, y),
                xytext=(dx, dy),
                textcoords="offset points",
                ha="right" if dx < 0 else "left",
                fontsize=9,
            )
        ax.set_title(titles[label], loc="left", fontsize=12)
        ax.set_xlabel("MDS 1 (inter-manifold distance)")
        ax.set_ylabel("MDS 2 (inter-manifold distance)")
        ax.text(0.01, 0.02, f"N={result['n_materials']}; dots=80% subsamples", transform=ax.transAxes, fontsize=8)
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(alpha=0.15, linewidth=0.5)
    fig.suptitle("Manifold-of-manifolds map: closer clouds preserve more similar material distances", fontsize=14)
    fig.savefig(ROOT / "figures/feature_manifold_distance_map.png", dpi=240)
    fig.savefig(ROOT / "figures/feature_manifold_distance_map.pdf")
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(SEED)
    main_jids, main_matrices = load_main_views()
    kl_jids, kl_matrices = load_kl_views()
    results = {
        "2d_atlas": bootstrap_map(main_matrices, "Structure", ["Transport n", "Transport p"], rng),
        "kl_validation": bootstrap_map(kl_matrices, "Electronic", ["Lattice kL"], rng),
    }
    points, pairs = [], []
    for panel, result in results.items():
        p, q = rows_for_export(panel, result)
        points.extend(p)
        pairs.extend(q)
    pd.DataFrame(pairs).to_csv(ROOT / "data/audit/manifold_of_manifolds_distances.csv", index=False)
    payload = {
        "method": "rank-distance correlation; d=sqrt(2*(1-rho)); classical MDS; 300x 80% subsampling",
        "panels": {
            "2d_atlas": {"n_materials": len(main_jids), "jids": main_jids},
            "kl_validation": {"n_materials": len(kl_jids), "jids": kl_jids},
        },
        "points": points,
        "pairs": pairs,
    }
    (ROOT / "data/processed/manifold_of_manifolds_map.json").write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    draw_static(results)
    print(pd.DataFrame(pairs).sort_values(["panel", "distance"]).to_string(index=False))
    print(f"\n2D atlas common N={len(main_jids)}; kL common N={len(kl_jids)}")
    print("Wrote PNG/PDF, distance table, and visualization JSON.")


if __name__ == "__main__":
    main()
