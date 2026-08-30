"""生成审计后的关键结果汇总图与机器可读表。"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config


def main():
    targets = pd.read_parquet(config.PROC_DIR / "kappa_L_targets.parquet")
    exp = targets[(targets["method"] == "experimental") & targets["material_id"].notna()]
    exp = exp.set_index("material_id")["kappa_L"]
    model = targets[targets["method"] == "Snyder-300K-model"].set_index("material_id")["kappa_L"]
    common = exp.index.intersection(model.index)
    x = np.log10(model.loc[common].to_numpy(float))
    y = np.log10(exp.loc[common].to_numpy(float))
    rho = float(stats.spearmanr(x, y).statistic)

    inc = pd.read_csv(config.PROC_DIR / "geometry_increment_folds.csv")
    inc_sum = pd.read_csv(config.PROC_DIR / "geometry_increment_summary.csv")
    overlap = pd.read_csv(config.PROC_DIR / "view_overlap.csv").set_index("pair")
    eg_pairs = ["Eg vs kL_clarke", "Eg vs kL_snyder"]

    rows = [
        {"result": "unique_experimental_mpid_matches", "value": len(common), "detail": "formula maps to one MP id"},
        {"result": "experimental_vs_snyder_spearman", "value": rho, "detail": f"N={len(common)}"},
    ]
    for _, r in inc_sum.iterrows():
        rows.append({"result": f"geometry_delta_r2_{r.target}", "value": r.delta_mean,
                     "detail": f"fold range [{r.delta_min:.4f}, {r.delta_max:.4f}], N={int(r.N)}"})
    for pair in eg_pairs:
        r = overlap.loc[pair]
        rows.append({"result": pair.replace(" ", "_"), "value": r["overlap"],
                     "detail": f"null={r.null_mean:.6f}, enrichment={r.enrichment:.3f}, z={r.z:.2f}"})
    pd.DataFrame(rows).to_csv(config.PROC_DIR / "corrected_key_results.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    ax.scatter(x, y, s=28, alpha=0.7, color="#2a6fbb", edgecolor="none")
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], "--", color="0.35", lw=1)
    ax.set_xlabel("log10 Snyder 300 K model κL")
    ax.set_ylabel("log10 experimental κL")
    ax.set_title(f"A  Unique formula→MP matches: N={len(common)}, ρ={rho:+.2f}")

    ax = axes[0, 1]
    order = ["snyder_acoustic", "clarke", "kappa_exp"]
    colors = {"snyder_acoustic": "#8d99ae", "clarke": "#e9c46a", "kappa_exp": "#2a9d8f"}
    for i, target in enumerate(order):
        vals = inc.loc[inc["target"] == target, "delta_r2_geometry"].to_numpy()
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.08, 0.08, len(vals)), vals,
                   color=colors[target], s=45, zorder=3)
        ax.plot([i - 0.18, i + 0.18], [vals.mean(), vals.mean()], color="black", lw=2)
    ax.axhline(0, color="0.3", lw=1)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Snyder model", "Clarke min", "Experiment\n(N=59)"])
    ax.set_ylabel("paired ΔR² from adding geometry SOAP")
    ax.set_title("B  Geometry increment: experiment crosses zero")

    ax = axes[1, 0]
    labels = ["Eg–Clarke", "Eg–Snyder model"]
    observed = np.array([overlap.loc[p, "overlap"] for p in eg_pairs])
    null = np.array([overlap.loc[p, "null_mean"] for p in eg_pairs])
    ci_lo = np.array([overlap.loc[p, "ci_lo"] for p in eg_pairs])
    ci_hi = np.array([overlap.loc[p, "ci_hi"] for p in eg_pairs])
    xx = np.arange(2)
    ax.bar(xx - 0.18, observed, width=0.36, label="observed", color="#457b9d")
    ax.bar(xx + 0.18, null, width=0.36, label="null", color="#adb5bd")
    ax.errorbar(xx - 0.18, observed, yerr=[observed - ci_lo, ci_hi - observed],
                fmt="none", ecolor="black", capsize=3)
    ax.set_xticks(xx); ax.set_xticklabels(labels)
    ax.set_ylabel("kNN overlap (k=10)")
    ax.set_title("C  kNN recomputed inside N=2472 common cohort")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.axis("off")
    statuses = [
        ("Q1  Geometry → true κL", "INCONCLUSIVE", "N=59; ΔR² fold range crosses 0", "#e9c46a"),
        ("Q2  κe ↔ true κL relation", "NOT IDENTIFIED", "N=70 formula join; 95% CI crosses zero", "#e76f51"),
        ("Q3  candidate ranking", "BLOCKED", "missing same-T by-id data, κe/τ, stability", "#d1495b"),
    ]
    for i, (question, status, why, color) in enumerate(statuses):
        y0 = 0.83 - i * 0.30
        ax.add_patch(plt.Rectangle((0.02, y0 - 0.10), 0.96, 0.22, color=color, alpha=0.18,
                                   transform=ax.transAxes))
        ax.text(0.05, y0 + 0.04, question, transform=ax.transAxes, fontsize=11, weight="bold")
        ax.text(0.95, y0 + 0.04, status, transform=ax.transAxes, ha="right",
                fontsize=11, weight="bold", color=color)
        ax.text(0.05, y0 - 0.045, why, transform=ax.transAxes, fontsize=10, color="0.25")
    ax.set_title("D  What the current data can support", loc="left")

    fig.suptitle("kappaL-refactored — corrected audit summary", fontsize=15, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    path = config.FIG_DIR / "corrected_audit_summary.png"
    plt.savefig(path, dpi=180)
    print(f"saved {path}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
