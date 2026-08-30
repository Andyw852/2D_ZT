#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""入口：运行多智能体材料发现闭环(默认 dry_run)。

用法:
  python run_pipeline.py                # 跑满 config.round_limit 轮
  python run_pipeline.py --rounds 1     # 只跑 1 轮
  python run_pipeline.py --submit       # 关闭 dry_run，真实提交(慎用)
"""
import os, sys, argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matexplore.config import load_config
from matexplore.orchestrator import run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=None)
    ap.add_argument("--submit", action="store_true", help="关闭 dry_run，真实提交到超算")
    args = ap.parse_args()
    cfg = load_config()
    if args.rounds:
        cfg["project"]["round_limit"] = args.rounds
    if args.submit:
        cfg["validation"]["dry_run"] = False
        print("WARNING: dry_run 已关闭，将真实提交作业到超算。", file=sys.stderr)
    run_id, results = run(cfg)
    print("RUN_ID:", run_id)
    print("results written under", os.path.join(cfg.store.root, run_id))


if __name__ == "__main__":
    main()

