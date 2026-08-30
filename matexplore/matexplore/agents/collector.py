# -*- coding: utf-8 -*-
"""智能体 4 —— 收集结果 + 写回 ledger(胜者/验证状态) + 反馈给下一轮。"""
import os

from .base import Agent


class CollectorAgent(Agent):
    name = "collector"

    def run(self, ctx):
        ledger = ctx.ledger
        cands = ctx.artifacts.get("validation", [])
        out = ctx.round_dir()

        # 本轮胜者：T0 验证后按 merit 排序取前若干，作为下一轮替换种子
        ranked = sorted(cands, key=lambda x: -x.get("merit", 0))
        winners = []
        for c in ranked[:10]:
            winners.append({
                "formula": c["formula"],
                "merit": c.get("merit"),
                "prediction": c.get("validated_prediction", {}),
                "species": c.get("species"),
                "lattice": c.get("lattice") if "lattice" in c else None,
                "coords": c.get("coords") if "coords" in c else None,
                "seed": c.get("seed_formula"),
            })
        ledger.set_winners(winners)

        # 记录验证状态(T0 已筛选)
        for c in ranked:
            ledger.mark_validated(c["formula"], "T0_screened",
                                  {"merit": c.get("merit"),
                                   "p_ZT_e": c.get("validated_prediction", {}).get("p_ZT_e")})

        top_summary = []
        for c in ranked[:10]:
            p = c["validated_prediction"]
            top_summary.append({
                "formula": c["formula"],
                "p_ZT_e": round(p.get("p_ZT_e", 0), 2),
                "n_ZT_e": round(p.get("n_ZT_e", 0), 2),
                "p_logPF": round(p.get("p_log10_PF", 0), 2),
                "seed": c.get("seed_formula"),
            })

        feedback = {
            "agent": self.name,
            "round": ctx.round,
            "n_validated": len(ranked),
            "winners_forwarded_to_next_round": [w["formula"] for w in winners],
            "top_candidates": top_summary,
            "next_round_actions": [
                "以本轮胜者为种子继续元素族替换(生成器已接 ledger)",
                "对最终胜者跑 opt-mace-gpu@3090 确认稳定性",
                "对稳定者跑 band-dft-cpu+ke-dft-cpu@jzzn 取 Eg/m*/S/σ/PF 真值",
                "用真值回填重算相关，更新假设权重(闭环)",
            ],
        }
        ctx.store.write(ctx.round, "feedback", feedback)
        ctx.store.write(ctx.round, "summary", {
            "round": ctx.round,
            "top_candidates": top_summary,
            "winners": [w["formula"] for w in winners],
            "feedback": feedback,
        })
        return feedback

