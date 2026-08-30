# -*- coding: utf-8 -*-
"""智能体 1 —— 分析/假设提出（机制优先，ARPS P1 语法）。

读既有数据库 -> 计算"易算特征<->热电目标"相关面板(证据) -> 生成机制优先假设：
  每条假设 = 设计变量 -> 可算描述子链 -> 目标性质 -> 方向 (+ 边界 + 可证伪条件 + next_best_test)。
纯统计相关只作为 evidence，不作为 claim 本身。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .base import Agent
from ..knowledge import database
from ..hypothesis.correlation import compute_correlation, plot_polyline


def _rho(corr, feature, target, carrier):
    r = corr[(corr.feature == feature) & (corr.target == target) & (corr.carrier == carrier)]
    return float(r.spearman.iloc[0]) if len(r) else None


def form_mechanisms(corr, panel):
    """从相关面板 + 物理知识生成机制优先假设。"""
    chi_eg_n = spearmanr(*[panel[k].dropna() for k in ("electronegativity_mean", "Eg_optb88vdw")]).statistic if False else None
    # 组成特征 <-> Eg 的混淆度(证明"组成特征是 Eg 的代理")
    def rho2(a, b):
        d = panel[[a, b]].dropna()
        return float(spearmanr(d[a], d[b]).statistic) if len(d) > 30 else None
    chi_eg = rho2("electronegativity_mean", "Eg_optb88vdw")
    ie_eg = rho2("ionization_energy_mean", "Eg_optb88vdw")

    M = [
        dict(id="H1", variable="m* (有效质量)",
             chain="m* → Seebeck S (Pisarenko/Mott) → ZT_e=S²/L",
             target="ZT_e", direction="+ (m*↑ ⇒ ZT_e↑)",
             evidence=dict(n=_rho(corr, "m* (median, m_e)", "ZT_e", "n"),
                           p=_rho(corr, "m* (median, m_e)", "ZT_e", "p")),
             boundary="T=600K, n=1e20, 二维单层, 常数弛豫时间",
             falsifier="同边界内某材料 m* 高但 ZT_e 低",
             next_best_test="固定骨架单元素替换，算 m*/Eg，看 S、ZT_e 是否同向变化",
             grade="screening(相关)→待 discrimination(受控对比)"),
        dict(id="H2", variable="Eg (带隙)",
             chain="Eg → 抑制双极输运 → S↑ → ZT_e↑",
             target="ZT_e", direction="+ (Eg↑ ⇒ ZT_e↑)",
             evidence=dict(n=_rho(corr, "Eg (band gap, eV)", "ZT_e", "n"),
                           p=_rho(corr, "Eg (band gap, eV)", "ZT_e", "p")),
             boundary="同上",
             falsifier="同边界内宽 Eg 但 ZT_e 低(应来自 m* 或弛豫时间异常)",
             next_best_test="对 seed 与替换后候选算 Eg，看方向",
             grade="screening→discrimination"),
        dict(id="H3", variable="电负性均值(组成)",
             chain="电负性 → 键共价性/带隙 → S → ZT_e (组成是 Eg 的代理)",
             target="ZT_e", direction="+ (电负性↑ ⇒ ZT_e↑)",
             evidence=dict(n=_rho(corr, "electronegativity mean", "ZT_e", "n"),
                           p=_rho(corr, "electronegativity mean", "ZT_e", "p"),
                           chi_vs_Eg=chi_eg, ie_vs_Eg=ie_eg),
             boundary="组成特征仅作无-DFT 粗筛",
             falsifier="若控制 Eg/m* 后组成特征无独立增量，则是纯代理",
             next_best_test="控制 Eg 的偏相关；或对同 Eg 不同组成算 S",
             grade="screening(相关，可能被 Eg 混淆)"),
        dict(id="H4", variable="m* (S–σ 权衡)",
             chain="m*↑ ⇒ σ↓ (μ∝1/m*) 但 PF=S²σ 仍随 S 上升",
             target="log10 PF", direction="PF 峰值在中等 m*/Eg(非单调)",
             evidence=dict(m_vs_sigma_n=_rho(corr, "m* (median, m_e)", "log10 sigma", "n"),
                           zte_vs_pf_n=_rho(corr, "m* (median, m_e)", "log10 PF", "n")),
             boundary="固定 n=1e20",
             falsifier="若重质量材料 PF 反而更高且 σ 也高，则权衡不成立",
             next_best_test="seed vs 替换后算 σ(PF)，看 S-σ 是否此消彼长",
             grade="screening→discrimination"),
        dict(id="H5", variable="层厚(真空/√面积)",
             chain="更薄二维层 → 量子限域 → S/ZT_e↑",
             target="ZT_e", direction="− (层厚↓ ⇒ ZT_e↑)",
             evidence=dict(n=_rho(corr, "vacuum / sqrt(area)", "ZT_e", "n"),
                           p=_rho(corr, "vacuum / sqrt(area)", "ZT_e", "p")),
             boundary="二维单层(真空轴足够大)",
             falsifier="若限域方向相反或层厚与 Eg 共变，则非独立机制",
             next_best_test="固定组成，改变真空/层厚算 S",
             grade="screening"),
    ]
    return M


class AnalystAgent(Agent):
    name = "analyst"

    def run(self, ctx):
        panel = database.load_panel(self.cfg)
        corr = compute_correlation(panel)
        mechanisms = form_mechanisms(corr, panel)

        out_dir = ctx.round_dir()
        os.makedirs(out_dir, exist_ok=True)
        corr.to_csv(os.path.join(out_dir, "correlation.csv"), index=False)
        png = plot_polyline(corr, os.path.join(out_dir, "correlation_polyline.png"))

        result = {
            "agent": self.name,
            "n_materials": int(len(panel)),
            "n_feature_target_pairs": int(len(corr)),
            "mechanisms": mechanisms,
            "figure": os.path.abspath(png),
        }
        ctx.store.write(ctx.round, "hypothesis", result)
        ctx.artifacts["correlation"] = corr
        ctx.artifacts["hypothesis"] = result
        return result

