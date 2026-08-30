"""Step 11：双通道数据可用性审计 + 公式级探索性分析。

zT = S²σT / (κ_e + κ_L)。好的热电材料需要 PF 高且 κ_L 低（phonon-glass electron-crystal）。
可检验版本需要同一 material_id、同温度的 PF 与 κ_L。当前本地数据不满足这些
条件，因此本脚本不再输出“重要性夹角=60°”这一伪精确结论。

数据（本地可用，诚实标注来源与温度）：
- 电子：JARVIS dft_3d 的 n/p-powerfact / n/p-Seebeck / ncond / nkappa（600 K、10²⁰ cm⁻³）
  —— 按 reduced_formula 映射到 MP（Ricci by-id 数据不在本地，见局限）。
- 晶格：MP snyder_acoustic（Snyder 300 K 解析模型，按 MP id）。

JARVIS↔MP 只能按 reduced_formula 连接，无法识别多晶型。本脚本先把 MP 多晶型
聚合到每个公式一行，再做带明确 ``formula_level_proxy`` 标签的探索性相关；它不能
回答材料级双通道是否可独立调控。
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config
from pymatgen.core import Composition

JARVIS = config.EXTERNAL_DATA_DIR / "jarvis_kl" / "jdft_3d-8-18-2021.json"


def _canon(f):
    try:
        return Composition(str(f)).reduced_formula
    except Exception:
        return None


def _finite(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return np.nan
    return np.nan if v <= -99998 else v


def build_table():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    meta["canon"] = meta["formula"].map(_canon)

    jar = json.load(open(JARVIS))
    best = {}
    for x in jar:
        canon = _canon(x.get("formula", ""))
        fe = _finite(x.get("formation_energy_peratom"))
        if not canon or np.isnan(fe):
            continue
        cur = best.get(canon)
        if cur is None or fe < cur["fe"]:
            best[canon] = {
                "fe": fe,
                "PF_n": _finite(x.get("n-powerfact")),
                "PF_p": _finite(x.get("p-powerfact")),
                "S_n": _finite(x.get("n-Seebeck")),
                "S_p": _finite(x.get("p-Seebeck")),
                "sigma_n": _finite(x.get("ncond")),
                "sigma_p": _finite(x.get("pcond")),
                "kappa_e_n": _finite(x.get("nkappa")),
                "kappa_e_p": _finite(x.get("pkappa")),
                "Eg_opt": _finite(x.get("optb88vdw_bandgap")),
            }

    jar_df = pd.DataFrame([{"canon": canon, **vals} for canon, vals in best.items()])
    # 一条公式对应多个 MP 多晶型时，保留中位数和 IQR，而不是复制同一个 PF 多次。
    agg = meta.groupby("canon", dropna=True).agg(
        formula=("formula", "first"),
        n_mp_polymorphs=("material_id", "size"),
        kappa_L_snyder_median=("snyder_acoustic", "median"),
        kappa_L_snyder_q25=("snyder_acoustic", lambda x: x.quantile(0.25)),
        kappa_L_snyder_q75=("snyder_acoustic", lambda x: x.quantile(0.75)),
        kappa_L_clarke_median=("clarke", "median"),
        bulk_vrh=("bulk_vrh", "median"),
        shear_vrh=("shear_vrh", "median"),
        debye=("debye", "median"),
        density=("density", "median"),
    ).reset_index()
    df = agg.merge(jar_df, on="canon", how="inner", validate="one_to_one")
    df["match_quality"] = "formula_level_proxy"
    df["T_electronic"] = 600.0     # JARVIS powerfact 参考条件
    df["T_kappa"] = 300.0
    df["doping"] = "1e20 cm-3"
    df.to_parquet(config.PROC_DIR / "dual_channel.parquet", index=False)
    print("dual-channel 公式级代理表（每个 canon 一行）:", len(df), "条")
    print("含多 MP 多晶型的公式占比:", round(float((df["n_mp_polymorphs"] > 1).mean()), 3))
    return df


def formula_proxy_association(df):
    """每个公式一票的探索性 Spearman + formula bootstrap CI。"""
    rng = np.random.RandomState(config.SEED)
    rows = []
    for carrier in ["n", "p"]:
        dfc = df[["canon", f"PF_{carrier}", "kappa_L_snyder_median"]].copy()
        dfc["logPF"] = np.log10(dfc[f"PF_{carrier}"])
        dfc["logkL"] = np.log10(dfc["kappa_L_snyder_median"])
        dfc = dfc.replace([np.inf, -np.inf], np.nan).dropna(subset=["logPF", "logkL"])
        rho = float(stats.spearmanr(dfc["logPF"], dfc["logkL"]).statistic)
        boots = []
        for _ in range(config.N_BOOTSTRAP):
            take = rng.randint(0, len(dfc), len(dfc))
            boots.append(stats.spearmanr(
                dfc["logPF"].to_numpy()[take], dfc["logkL"].to_numpy()[take]).statistic)
        lo, hi = np.quantile(boots, [0.025, 0.975])
        rows.append({
            "carrier": carrier,
            "N_unique_formula": len(dfc),
            "spearman_logPF_logkL": rho,
            "ci_lo": float(lo),
            "ci_hi": float(hi),
            "analysis_status": "exploratory_formula_level_only",
        })
    result = pd.DataFrame(rows)
    result.to_csv(config.PROC_DIR / "dual_channel_formula_corr.csv", index=False)
    print(result.to_string(index=False))

    # 公式级代理散点；不再绘制没有不确定性的 RF 重要性夹角。
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for ax, carrier in zip(axes, ["n", "p"]):
        d = df[[f"PF_{carrier}", "kappa_L_snyder_median", "n_mp_polymorphs"]].copy()
        d = d[(d[f"PF_{carrier}"] > 0) & (d["kappa_L_snyder_median"] > 0)]
        ax.scatter(np.log10(d["kappa_L_snyder_median"]), np.log10(d[f"PF_{carrier}"]),
                   s=8, alpha=0.25, c=np.clip(d["n_mp_polymorphs"], 1, 10), cmap="viridis")
        rr = result[result["carrier"] == carrier].iloc[0]
        ax.set_title(f"{carrier}-type: ρ={rr.spearman_logPF_logkL:+.2f} "
                     f"[{rr.ci_lo:+.2f}, {rr.ci_hi:+.2f}]")
        ax.set_xlabel("log10 Snyder κL model (MP formula median)")
        ax.set_ylabel("log10 JARVIS PF")
    fig.suptitle("Formula-level proxy only — not a material-id matched dual-channel test")
    plt.tight_layout(); plt.savefig(config.FIG_DIR / "dual_channel_formula_proxy.png", dpi=160)
    print("saved figures/dual_channel_formula_proxy.png")
    return result


def main():
    df = build_table()
    print()
    print("=== 公式级探索性关联（不能据此判断双通道正交性）===")
    formula_proxy_association(df)


if __name__ == "__main__":
    main()
