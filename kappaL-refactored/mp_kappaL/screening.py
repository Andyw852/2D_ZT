"""Step 12：筛选前置条件验收。

只有同一 material_id、同温度、可比较单位的 PF/κe/κL 与稳定性数据齐备时，
才允许产生 Pareto 前沿或候选排名。当前本地数据不满足条件，脚本会输出可审计的
readiness 表和阻塞报告，而不是继续发布误导性的 ``PF/κL`` 排名。

旧的 ``PF / Snyder-κL`` 排序把公式级电子输运复制到 MP 多晶型，混合 600/300 K，
并忽略 κe、τ 和稳定性，已全部撤回；已知材料名次也不再当作验收证据。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config


def pareto_mask(pf, kl):
    pf = np.asarray(pf); kl = np.asarray(kl)
    keep = np.ones(len(pf), dtype=bool)
    for i in range(len(pf)):
        dominates = (pf >= pf[i]) & (kl <= kl[i]) & ((pf > pf[i]) | (kl < kl[i]))
        if dominates.any():
            keep[i] = False
    return keep


def main():
    df = pd.read_parquet(config.PROC_DIR / "dual_channel.parquet")
    checks = pd.DataFrame([
        ("electronic_and_kappa_same_material_id", False,
         "JARVIS↔MP only reduced_formula; polymorph identity is unknown"),
        ("same_temperature", False, "PF at 600 K; Snyder model fixed at 300 K"),
        ("independent_kappa_L_target", False,
         "large-N target is an elastic-property formula, not measured/BTE κL"),
        ("kappa_e_units_and_tau_validated", False,
         "local JARVIS fields are not converted into a validated zT denominator"),
        ("energy_above_hull_available", False, "stability field is absent locally"),
        ("unique_formula_rows", bool(df["canon"].is_unique),
         "formula proxy table is de-duplicated, but this does not resolve polymorph mapping"),
    ], columns=["criterion", "passed", "evidence"])
    checks.to_csv(config.PROC_DIR / "screening_readiness.csv", index=False)

    # 清空旧的伪候选结果，保留显式 schema 与阻塞原因，避免下游误读陈旧排名。
    pd.DataFrame(columns=["rank", "material_id", "formula", "status", "reason"]).to_csv(
        config.PROC_DIR / "candidates.csv", index=False)

    lines = [
        "# 热电候选筛选状态（Step 12）",
        "",
        "> **结论：当前数据不允许生成候选排名。** 旧版 `PF_best / κ_L_snyder` 排名已撤回。",
        "",
        "## 前置条件验收",
        "",
        "| 条件 | 通过 | 证据 |",
        "|---|---:|---|",
    ]
    for _, r in checks.iterrows():
        lines.append(f"| {r['criterion']} | {'是' if r['passed'] else '否'} | {r['evidence']} |")
    lines += [
        "",
        "## 为什么撤回旧排名",
        "",
        "旧排名把 600 K 的 JARVIS PF 与 300 K 的 Snyder 解析模型按化学式拼接，",
        "把同一 PF 复制给多个 MP 多晶型，忽略 κe、τ、energy_above_hull，再把 PF/κL",
        "称作 zT 排序代理。该顺序不能由 zT 公式推出，已知材料的名次也不能作为验证。",
        "",
        "## 恢复筛选所需数据",
        "",
        "- material_id 对齐的电子输运（PF、σ/τ、κe/τ）与载流子浓度、温度；",
        "- 同一结构、同一温度的实验或 BTE κL；",
        "- energy_above_hull；",
        "- τ 敏感性和已知热电体系的盲法召回验证。",
    ]
    (config.REPORTS_DIR / "candidates.md").write_text("\n".join(lines), encoding="utf-8")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.8))
    y = np.arange(len(checks))
    colors = np.where(checks["passed"], "#2a9d8f", "#d1495b")
    ax.barh(y, np.ones(len(checks)), color=colors)
    ax.set_yticks(y); ax.set_yticklabels(checks["criterion"])
    ax.set_xlim(0, 1); ax.set_xticks([])
    ax.invert_yaxis()
    for i, passed in enumerate(checks["passed"]):
        ax.text(0.5, i, "PASS" if passed else "BLOCKED", ha="center", va="center",
                color="white", fontweight="bold")
    ax.set_title("Screening readiness: candidate ranking is blocked")
    plt.tight_layout(); plt.savefig(config.FIG_DIR / "screening_readiness.png", dpi=160)
    print(checks.to_string(index=False))
    print("saved screening_readiness.csv, empty candidates.csv, candidates.md, screening_readiness.png")


if __name__ == "__main__":
    main()
