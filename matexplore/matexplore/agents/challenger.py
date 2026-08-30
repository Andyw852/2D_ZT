# -*- coding: utf-8 -*-
"""智能体 1.5 —— 对抗性审稿(Challenger，ARPS 借鉴)。

独立批评分析师提出的机制假设：攻其边界/混淆/替代机制/决策价值，
并给出下一轮该跑的甄别性实验。不执行计算、不写结论，只出 issues + 推荐。
"""
import os
import numpy as np
from scipy.stats import spearmanr

from .base import Agent
from ..knowledge import database


def _rho(panel, a, b):
    d = panel[[a, b]].dropna()
    if len(d) < 30:
        return None
    return float(spearmanr(d[a], d[b]).statistic)


class ChallengerAgent(Agent):
    name = "challenger"

    def run(self, ctx):
        hyp = ctx.artifacts.get("hypothesis", {})
        mechanisms = hyp.get("mechanisms", [])
        panel = database.load_panel(self.cfg)

        chi_eg = _rho(panel, "electronegativity_mean", "Eg_optb88vdw")
        ie_eg = _rho(panel, "ionization_energy_mean", "Eg_optb88vdw")
        area_eg = _rho(panel, "inplane_area", "Eg_optb88vdw")
        chi_m = _rho(panel, "electronegativity_mean", "m_elec_median")

        issues = []
        for m in mechanisms:
            mid = m["id"]
            needs_dft = ("m*" in m["variable"]) or ("Eg" in m["variable"])
            alt = None
            if mid == "H3":
                alt = ("组成特征(电负性/电离能)可能是 Eg 的代理而非独立机制: chi_Eg=%.2f, IE_Eg=%.2f") % (chi_eg or 0, ie_eg or 0)
            elif mid == "H5":
                alt = "层厚与 Eg/组成可能共变，非独立限域机制"
            elif mid in ("H1", "H2"):
                alt = "m* 与 Eg 强相关(+0.7+)，二者不可分：H1/H2 可能是同一机制的两种描述"
            single_var = ("单元素替换同时改变电负性/原子质量/原子半径(如 Pb换Sn)，非严格单变量，需声明为多描述子联动对比")
            if needs_dft:
                issues.append(dict(
                    target_id=mid, severity="high", issue_type="methodology",
                    text="%s 需要 DFT 实测；当前只有相关证据(screening 级)，不能当机制结论" % m["variable"],
                    consequence="若不加甄别性实验，会把相关误当因果，推荐候选可能是假阳性",
                    location="within-boundary",
                ))
            else:
                issues.append(dict(
                    target_id=mid, severity="low", issue_type="boundary",
                    text="纯组成/几何特征的可算链成立，但需声明代理关系",
                    consequence="若混淆未排除，会把组成代理误当独立机制",
                    location="within-boundary",
                ))
            if alt:
                issues.append(dict(
                    target_id=mid, severity="high" if mid == "H3" else "medium",
                    issue_type="statistical", text="最强未排除替代机制: " + alt,
                    consequence="若成立，" + mid + " 的方向可能是伪的",
                    location="within-boundary",
                ))
            issues.append(dict(
                target_id=mid, severity="medium", issue_type="methodology",
                text=single_var,
                consequence="替换前后对比是联动对比而非受控单变量，结论需限定在给定骨架内",
                location="within-boundary",
            ))

        recommended = {
            "abstract": ("对 top 候选与其 seed 做受控对比：算 Eg(及 m*)，看替换是否按预测方向移动 Eg，从而甄别 H2/H3 是因果还是代理"),
            "feasible_on_3090": [
                "ML 带隙预测(装 MatGL，pip 可装)：对 seed 与替换后候选预测 Eg，看方向",
                "声子动力学稳定性(phonon-mace，MACE+phonopy 已在 3090)：查虚频，过滤动力学不稳定候选(任何电子性质讨论之前的第一道真值闸)",
                "更强甄别需 DFT(OpenMX/QE 可装)：PBE 带隙 + 有效质量",
            ],
            "decision_value": "把相关暗示升级为受控对比结论，直接改变候选排名",
        }

        result = {
            "agent": self.name,
            "confounding": {"chi_vs_Eg": chi_eg, "ie_vs_Eg": ie_eg,
                            "area_vs_Eg": area_eg, "chi_vs_m": chi_m},
            "issues": issues,
            "readiness": "needs-discrimination",
            "recommended_discriminating_test": recommended,
            "recommended_next_round": ("跑一轮甄别：3090 上对 top 候选做 (a)声子动力学稳定性 + (b)ML 带隙 的 seed-vs-替换 对比"),
        }
        ctx.store.write(ctx.round, "challenge", result)
        ctx.artifacts["challenge"] = result
        return result
