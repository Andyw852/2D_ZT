# -*- coding: utf-8 -*-
"""智能体 3 —— 验证生成材料的性质。

三层验证(由浅入深)：
  T0 替代筛选(te-screen)：易算特征 -> ZT_e/PF 预测，登录节点秒级，任何机器可跑。
  T1 MACE 弛豫(opt-mace-gpu @3090)：稳定性/形成能，GPU 快速。
  T2 DFT 电子输运(band-dft-cpu + ke-dft-cpu @jzzn)：Eg/m*/S/σ/PF/κe 真值。
(晶格热导率 kl-* 本轮跳过。)

本实现默认 dry_run=true：本地完成 T0 独立复算 + 写出 POSCAR + 生成 T1/T2 提交计划，
不真正 sbatch；置 dry_run=false 即真实提交。
"""
import os

from .base import Agent
from ..generation.surrogate import SurrogateScorer
from ..hypothesis.cheap_features import load_element_properties
from ..validation.taskflow_client import TaskflowClient


class ValidatorAgent(Agent):
    name = "validator"

    def run(self, ctx):
        cfg = self.cfg
        candidates = ctx.artifacts.get("candidates", [])
        scorer = SurrogateScorer(cfg.knowledge.surrogate_models_dir)
        elem_props = load_element_properties(cfg.knowledge.element_properties)
        tc = TaskflowClient(cfg)

        out_dir = ctx.round_dir()
        poscar_dir = os.path.join(out_dir, "poscars")
        os.makedirs(poscar_dir, exist_ok=True)

        validated = []
        plans = []
        for c in candidates:
            # T0 独立复算(与生成器共用 cheap_features 但独立走一遍)
            screen = tc.screen_locally(c, elem_props, scorer)
            c["validated_prediction"] = screen["prediction"]
            c["validated_features"] = screen["features"]
            # 写 POSCAR
            fn = os.path.join(poscar_dir, f"{c['formula']}.vasp")
            with open(fn, "w") as f:
                f.write(c["poscar"])
            c["poscar_path"] = fn
            validated.append(c)
            # T1/T2 提交计划
            plans.append(tc.plan_submission(c, cfg.validation.skill_screen, "3090"))
            plans.append(tc.plan_submission(c, cfg.validation.skill_relax, "3090"))
            plans.append(tc.plan_submission(c, cfg.validation.skill_band, "jzzn"))
            plans.append(tc.plan_submission(c, cfg.validation.skill_transport, "jzzn"))

        result = {
            "agent": self.name,
            "dry_run": cfg.validation.dry_run,
            "n_validated": len(validated),
            "candidates": validated,
            "submission_plans": plans,
        }
        ctx.store.write(ctx.round, "validation", result)
        ctx.artifacts["validation"] = validated
        return result

