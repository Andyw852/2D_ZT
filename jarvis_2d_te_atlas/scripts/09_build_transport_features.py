"""Phase I: 构建 permutation-invariant 输运特征（n/p 分别）。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

def eigs(jid, field):
    if field not in raw.get(jid, {}):
        return None
    return np.array([float(x) for x in raw[jid][field]], dtype=float)

def seebeck_feats(e):
    if e is None:
        return {}
    m = e.mean(); s = e.std(); mn = e.min(); mx = e.max()
    abs_e = np.abs(e)
    return {
        "S_mean": m, "S_std": s, "S_min": mn, "S_max": mx, "S_range": mx - mn,
        "S_abs_mean": abs_e.mean(), "S_abs_max": abs_e.max(),
        "S_relative_spread": s / (abs(m) + 1e-12),
    }

def pos_feats(e, eps_log, prefix):
    """sigma / kappa_e / PF 的正值 permutation-invariant 特征。"""
    if e is None:
        return {}
    m = e.mean(); s = e.std(); mn = e.min(); mx = e.max()
    anis = np.log10(mx / mn) if mn > 1e-9 else np.nan
    out = {
        prefix + "_mean": m, prefix + "_std": s, prefix + "_min": mn, prefix + "_max": mx,
        prefix + "_anisotropy_log": anis,
    }
    return out

# epsilon for log transform (第 29 节: min_positive / 10)
all_sigma_mean = [np.mean([float(x) for x in raw[j]["ncond"]]) for j in raw if "ncond" in raw[j]]
all_sigma_mean += [np.mean([float(x) for x in raw[j]["pcond"]]) for j in raw if "pcond" in raw[j]]
all_kappa_mean = [np.mean([float(x) for x in raw[j]["nkappa"]]) for j in raw if "nkappa" in raw[j]]
all_kappa_mean += [np.mean([float(x) for x in raw[j]["pkappa"]]) for j in raw if "pkappa" in raw[j]]
eps_sigma = min(v for v in all_sigma_mean if v > 0) / 10
eps_kappa = min(v for v in all_kappa_mean if v > 0) / 10
print(f"eps_sigma = {eps_sigma:.6g}, eps_kappa = {eps_kappa:.6g}")

def build(carrier):
    n = "n" if carrier == "n" else "p"
    rows = []
    for jid in raw:
        eS = eigs(jid, n + "seeb")
        if eS is None and eigs(jid, n + "cond") is None:
            continue
        row = {"jid": jid}
        row.update(seebeck_feats(eS))
        row.update(pos_feats(eigs(jid, n + "cond"), eps_sigma, "sigma"))
        row.update(pos_feats(eigs(jid, n + "kappa"), eps_kappa, "kappa_e"))
        row.update(pos_feats(eigs(jid, n + "pf"), 0, "PF"))
        # log transforms (第 29 节)
        sm = row.get("sigma_mean"); km = row.get("kappa_e_mean")
        row["log_sigma_mean"] = np.log10(sm + eps_sigma) if sm is not None and not np.isnan(sm) else np.nan
        row["log_kappa_e_mean"] = np.log10(km + eps_kappa) if km is not None and not np.isnan(km) else np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    df = df.sort_values("jid").reset_index(drop=True)
    return df

for carrier in ["n", "p"]:
    df = build(carrier)
    out = root / "features" / "transport" / f"{carrier}_transport_tensor_features.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"{carrier}-type: {df.shape[0]} rows, {df.shape[1]} cols -> {out.name}")
    print("  columns:", list(df.columns))
    print("  non-null counts:", int(df[["S_mean","sigma_mean","kappa_e_mean","PF_mean"]].notna().all(axis=1).sum()))
