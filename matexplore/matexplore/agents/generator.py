# -*- coding: utf-8 -*-
"""智能体 2 —— 材料生成(跨轮有记忆)。

第 0 轮：以面板中 p 型 ZT_e 最高的若干材料为种子，做元素族同价替换。
第 N 轮：额外加入"上一轮胜者"的结构作种子(围绕胜者继续替换)，并跳过
ledger 里已生成过的公式(去重)，实现"围绕胜者迭代精炼"的闭环。
"""
import os
import numpy as np

from .base import Agent
from ..knowledge import database
from ..generation.structure_models import (SubstitutionGenerator, make_poscar,
                                           formula_from_species)
from ..generation.surrogate import SurrogateScorer
from ..hypothesis.cheap_features import load_element_properties, compute_features


class GeneratorAgent(Agent):
    name = "generator"

    def run(self, ctx):
        cfg = self.cfg
        snapshot = database.load_snapshot(cfg)
        panel = database.load_panel(cfg)
        scorer = SurrogateScorer(cfg.knowledge.surrogate_models_dir)
        elem_props = load_element_properties(cfg.knowledge.element_properties)
        ledger = ctx.ledger
        seen = ledger.seen_formulas()

        # ---- 基线种子：面板 p 型 ZT_e 最高的材料结构 ---------------------
        top = (panel.sort_values("ZT_e_p", ascending=False)
               .dropna(subset=["ZT_e_p"]).head(60)["jid"].tolist())
        by_jid = {d["attributes"]["_jarvis_jid"]: d["attributes"] for d in snapshot}
        seed_structs = []
        for jid in top:
            a = by_jid.get(jid)
            if a is not None:
                seed_structs.append(a)
            if len(seed_structs) >= 30:
                break

        # ---- 上一轮胜者作为额外种子(闭环精炼) ---------------------------
        for w in ctx.artifacts.get("previous_winners", []):
            if w.get("species") and w.get("lattice") and w.get("coords"):
                seed_structs.append(dict(lattice_vectors=w["lattice"],
                                         species_at_sites=w["species"],
                                         cartesian_site_positions=w["coords"],
                                         _jarvis_formula=w.get("formula", "winner")))

        subg = SubstitutionGenerator(cfg.generation.substitution_groups)
        candidates = {}
        for a in seed_structs:
            base_sp = a["species_at_sites"]
            lv = a["lattice_vectors"]
            coords = a["cartesian_site_positions"]
            for new_sp in subg.expand(base_sp):
                formula = formula_from_species(new_sp)
                if formula in seen or formula in candidates:
                    continue
                poscar = make_poscar(lv, new_sp, coords, comment=formula)
                feat, meta = compute_features(poscar, elem_props)
                pred = scorer.predict(feat)
                candidates[formula] = dict(formula=formula, poscar=poscar,
                                           species=new_sp, prediction=pred,
                                           lattice=lv, coords=coords,
                                           seed_formula=a.get("_jarvis_formula"),
                                           seed_jid=a.get("_jarvis_jid", "winner"))

        # 评分排序取 top_k
        for c in candidates.values():
            c["merit_p"] = scorer.merit(c["prediction"], "p")
            c["merit_n"] = scorer.merit(c["prediction"], "n")
            c["merit"] = max(c["merit_p"], c["merit_n"])
        ranked = sorted(candidates.values(), key=lambda c: -c["merit"])
        topk = ranked[:cfg.project.top_k_per_round]

        # 记入 ledger(去重记忆)
        ledger.mark_seen(list(candidates.keys()))

        result = {
            "agent": self.name,
            "n_baseline_seeds": len(seed_structs) - len(ctx.artifacts.get("previous_winners", [])),
            "n_winner_seeds": len(ctx.artifacts.get("previous_winners", [])),
            "n_new_generated": len(candidates),
            "n_seen_skipped": 0,
            "top_k": topk,
        }
        ctx.store.write(ctx.round, "candidates", result)
        ctx.artifacts["candidates"] = topk
        return result

