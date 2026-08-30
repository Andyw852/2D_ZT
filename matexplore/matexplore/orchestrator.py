# -*- coding: utf-8 -*-
"""编排器：多智能体 假设->对抗审稿->生成->验证->反馈 闭环(ARPS 借鉴)。

每轮：analyst(机制假设) -> challenger(对抗审稿+推荐甄别实验)
          -> generator(受控替换生成) -> validator(分级验证) -> collector(汇总反馈)。
跨轮记忆走 Ledger(去重 + 胜者种子 + 验证状态 + 前沿)。
"""
import os
import sys

from .config import load_config
from .knowledge.round_store import RoundStore, Ledger
from .agents.analyst import AnalystAgent
from .agents.challenger import ChallengerAgent
from .agents.generator import GeneratorAgent
from .agents.validator import ValidatorAgent
from .agents.collector import CollectorAgent


class Ctx:
    def __init__(self, store, cfg, round_no, ledger):
        self.store = store
        self.cfg = cfg
        self.round = round_no
        self.ledger = ledger
        self.artifacts = {}

    def round_dir(self):
        return self.store.round_dir(self.round)


def build_agents(cfg):
    return [
        AnalystAgent(cfg),
        ChallengerAgent(cfg),
        GeneratorAgent(cfg),
        ValidatorAgent(cfg),
        CollectorAgent(cfg),
    ]


def run(cfg, run_id=None):
    store = RoundStore(cfg.store.root)
    if run_id is None:
        run_id = store.new_run()
    else:
        store.run_id = run_id
        store.run_dir = os.path.join(store.root, run_id)
    ledger = Ledger(store.run_dir)
    limit = cfg.project.round_limit
    results = []
    for rnd in range(limit):
        ctx = Ctx(store, cfg, rnd, ledger)
        ctx.artifacts["previous_winners"] = ledger.winners()
        round_result = {"round": rnd}
        for agent in build_agents(cfg):
            print("[round %d] %s ..." % (rnd, agent.name), file=sys.stderr)
            out = agent.run(ctx)
            round_result[agent.name] = (out or {})
        results.append(round_result)
    return run_id, results


def main():
    cfg = load_config()
    run_id, results = run(cfg)
    print("RUN_ID:", run_id)


if __name__ == "__main__":
    main()
