# -*- coding: utf-8 -*-
"""taskflow (tf) 客户端封装：技能发现 / 项目脚手架 / 提交命令生成。

验证器通过本模块衔接 taskflow：
  - 快速替代筛选  : tf -tt te-screen   -p <材料> start   (登录节点 run:gen)
  - MACE 弛豫(3090): tf -tt opt-mace-gpu -p <材料> start  (GPU)
  - DFT 电子输运(jzzn, 可选后续) : band-dft-cpu + ke-dft-cpu
晶格热导率(kl-*)本轮按用户要求跳过。
"""
import os
import subprocess


class TaskflowClient:
    def __init__(self, cfg):
        self.cfg = cfg
        self.tf = cfg.validation.tf_bin
        self.root = cfg.validation.taskflow_root
        self.dry_run = bool(cfg.validation.dry_run)

    def _run(self, args, cwd=None):
        cmd = [self.tf] + args
        if self.dry_run:
            return 0, "[dry-run] " + " ".join(cmd)
        p = subprocess.run(cmd, cwd=cwd or self.root, capture_output=True, text=True)
        return p.returncode, (p.stdout or "") + (p.stderr or "")

    def list_skills(self):
        code, out = self._run(["skills"])
        return code, out

    def plan_submission(self, candidate, skill, hpc=None):
        """为候选材料生成 taskflow 提交计划(命令 + 说明)。"""
        hpc = hpc or self.cfg.validation.default_hpc
        mat = candidate["formula"]
        lines = []
        lines.append(f"# 材料 {mat} —— 技能 {skill} 提交到 {hpc}")
        lines.append(f"tf -tt {skill} -p {mat} hpc {hpc}")
        lines.append(f"tf -tt {skill} -p {mat} -j 1 init   # 只生成输入不提交(先检查)")
        lines.append(f"tf -tt {skill} -p {mat} start      # 提交")
        return "\n".join(lines)

    def screen_locally(self, candidate, elem_props, scorer):
        """用 te-screen 技能同一套逻辑(纯 numpy)做替代打分，等价于 tf 技能产出。"""
        from ..hypothesis.cheap_features import compute_features
        feat, meta = compute_features(candidate["poscar"], elem_props)
        pred = scorer.predict(feat)
        return dict(features=feat, prediction=pred, meta=meta)

