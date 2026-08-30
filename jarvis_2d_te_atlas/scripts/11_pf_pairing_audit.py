"""L0-A: PF 本征值配对歧义审计。

PF = S_i^2 * sigma_p(i) / 1e6。由于 np.linalg.eigvals 分别作用于 S 与 sigma tensor，
S 与 sigma 的本征值顺序未必对应同一物理方向。本脚本枚举 sigma 的全部 6 种 permutation，
评估 PF 对配对顺序的敏感度。
"""
import json
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

def audit(carrier):
    n = "n" if carrier == "n" else "p"
    sf, cf, pf = n + "seeb", n + "cond", n + "pf"
    rows = []
    for jid in raw:
        if not all(k in raw[jid] for k in (sf, cf, pf)):
            continue
        S = np.array([float(x) for x in raw[jid][sf]])
        C = np.array([float(x) for x in raw[jid][cf]])
        PF_jarvis = np.array([float(x) for x in raw[jid][pf]])
        S2 = S ** 2
        pf_means = []
        for perm in itertools.permutations(range(3)):
            pf_i = S2 * C[list(perm)] / 1e6
            pf_means.append(pf_i.mean())
        pf_means = np.array(pf_means)
        mn, mx, med = pf_means.min(), pf_means.max(), np.median(pf_means)
        jarvis_mean = PF_jarvis.mean()
        # JARVIS pairing 对应 identity permutation
        identity = (S2 * C / 1e6).mean()
        rows.append({
            "jid": jid,
            "PF_jarvis_mean": jarvis_mean,
            "PF_identity_mean": identity,
            "PF_perm_min": mn, "PF_perm_max": mx, "PF_perm_median": med,
            "PF_perm_std": pf_means.std(),
        })
    df = pd.DataFrame(rows)
    # epsilon = min non-zero |median| / 10
    pos = df["PF_perm_median"].abs()
    pos = pos[pos > 0]
    eps = pos.min() / 10 if len(pos) else 1e-12
    df["ambiguity"] = (df["PF_perm_max"] - df["PF_perm_min"]) / (df["PF_perm_median"].abs() + eps)
    # JARVIS 相对位置：0=min, 1=max, 0.5=median
    df["jarvis_rel_pos"] = (df["PF_jarvis_mean"] - df["PF_perm_min"]) / (df["PF_perm_max"] - df["PF_perm_min"] + 1e-12)
    return df, eps

for carrier in ["n", "p"]:
    df, eps = audit(carrier)
    out = root / "data" / "audit" / f"pf_pairing_ambiguity_{carrier}.csv"
    df.to_csv(out, index=False)
    amb = df["ambiguity"]
    print(f"\n=== {carrier}-type PF pairing ambiguity (n={len(df)}, eps={eps:.3g}) ===")
    print(f"  median={amb.median():.4f} mean={amb.mean():.4f} p90={amb.quantile(0.9):.4f} p95={amb.quantile(0.95):.4f} max={amb.max():.4f}")
    for lo, hi, lab in [(0, 0.05, "<0.05"), (0.05, 0.20, "0.05-0.20"), (0.20, 0.50, "0.20-0.50"), (0.50, np.inf, ">=0.50")]:
        frac = ((amb >= lo) & (amb < hi)).mean()
        print(f"  ambiguity {lab}: {frac:.3f}")
    # JARVIS pairing 相对位置分布
    print(f"  JARVIS rel_pos: median={df['jarvis_rel_pos'].median():.3f} (0=min,1=max,0.5=median)")
    # 图
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].hist(np.log10(amb.clip(lower=1e-6)), bins=50, color="#4477aa")
    axes[0].set_xlabel("log10(PF pairing ambiguity)")
    axes[0].set_ylabel("count")
    axes[0].set_title(f"{carrier}-type PF pairing ambiguity")
    axes[1].scatter(df["PF_perm_median"], df["ambiguity"], s=4, alpha=0.5)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("PF_perm_median")
    axes[1].set_ylabel("ambiguity")
    axes[1].set_title(f"{carrier}-type ambiguity vs PF")
    plt.tight_layout()
    for suffix, idx in [("hist", 0)]:
        pass
    fig.savefig(root / "figures" / f"pf_pairing_ambiguity_hist_{carrier}.png", dpi=130)
    plt.close(fig)
    print(f"  wrote pf_pairing_ambiguity_{carrier}.csv + figures")
