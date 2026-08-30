# -*- coding: utf-8 -*-
"""闭环结果的持久化存储：每轮一个目录 + 一个跨轮 Ledger(已生成/胜者/验证状态)。

目录结构:
  runs/<run_id>/
    ledger.json         (跨轮累积: 已生成公式 / 胜者结构 / 验证状态)
    round_000/
      hypothesis.json   (分析智能体产出)
      candidates.json   (生成智能体产出)
      validation.json   (验证智能体产出)
      feedback.json     (收集智能体产出 -> 反馈给下一轮)
    state.json
"""
import json
import os
import time
import uuid


class RoundStore:
    def __init__(self, root):
        self.root = root
        self.run_id = None
        self.run_dir = None
        os.makedirs(root, exist_ok=True)

    def new_run(self):
        self.run_id = time.strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]
        self.run_dir = os.path.join(self.root, self.run_id)
        os.makedirs(self.round_dir(0), exist_ok=True)
        self._write_state(dict(round=0, status="running"))
        return self.run_id

    def round_dir(self, rnd):
        return os.path.join(self.run_dir, f"round_{rnd:03d}")

    def write(self, rnd, name, obj):
        d = self.round_dir(rnd)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name + ".json"), "w") as f:
            json.dump(obj, f, ensure_ascii=False, indent=2, default=str)

    def read(self, rnd, name):
        p = os.path.join(self.round_dir(rnd), name + ".json")
        if not os.path.exists(p):
            return None
        with open(p) as f:
            return json.load(f)

    def _write_state(self, state):
        with open(os.path.join(self.run_dir, "state.json"), "w") as f:
            json.dump(state, f, ensure_ascii=False, indent=2)


class Ledger:
    """跨轮累积记忆：已生成公式(去重) / 各轮胜者(供下一轮做种子) / 验证状态。"""

    def __init__(self, run_dir):
        self.path = os.path.join(run_dir, "ledger.json")
        self.data = {"seen": [], "winners": [], "validated": {}}
        if os.path.exists(self.path):
            with open(self.path) as f:
                self.data = json.load(f)

    def seen_formulas(self):
        return set(self.data["seen"])

    def mark_seen(self, formulas):
        s = set(self.data["seen"])
        s.update(formulas)
        self.data["seen"] = sorted(s)
        self.save()

    def winners(self):
        return self.data["winners"]

    def set_winners(self, winners):
        self.data["winners"] = winners
        self.save()

    def mark_validated(self, formula, status, info):
        self.data["validated"][formula] = {"status": status, **info}
        self.save()

    def validated(self):
        return self.data["validated"]

    def save(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with open(self.path, "w") as f:
            json.dump(self.data, f, ensure_ascii=False, indent=2, default=str)

