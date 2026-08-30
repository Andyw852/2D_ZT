#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""顺序验证(jzzn)：一次一个材料，MACE 链 opt-mace-cpu->phonon-mace-cpu->kl-mace-cpu。
unihamgnn(带隙, 3090)已单独收集，不在此链内。
"""
import os, sys, time, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tf_validator as tv

PROJECT_ROOT = tv.PROJECT_ROOT
tf_run = tv.tf_run
collect = tv.collect

# 仅 jzzn MACE 链(按依赖顺序)
JZZN_SKILLS = ["opt-mace-cpu", "phonon-mace-cpu", "kl-mace-cpu"]

# 已知基准（验证 POSCAR 转换 + MACE 链正确性，不是「发现」）：WSe2 有文献 κ_L 可对照
BENCHMARK_MATS = ["WSe2"]
# 候选结构（真正要探索/排序的材料）
CANDIDATE_MATS = ["MoS2", "SnSe2", "GaSe", "InSe", "Bi2Te3",
                  "F4Sn", "F4Ge", "ClF3Pb", "BrF3Pb", "BrCl2FPb2", "Br2ClFPb2"]
# 顺序验证列表：先基准后候选
MATS = BENCHMARK_MATS + CANDIDATE_MATS


def wait_done(skill, mat, minutes=35):
    deadline = time.time() + minutes * 60
    while time.time() < deadline:
        tf_run("-tt", skill, "-p", mat, "fetch", timeout=200)
        if collect(skill, mat) is not None:
            return True
        time.sleep(120)
    return False


def main():
    for mat in MATS:
        print("===== %s =====" % mat, flush=True)
        for skill in JZZN_SKILLS:
            d = collect(skill, mat)
            if d is not None:
                print("  [%s] 已有 %s" % (skill, json.dumps(d, ensure_ascii=False)[:140]), flush=True)
                continue
            for attempt in range(1, 3):
                rc, out, err = tf_run("-tt", skill, "-p", mat, "start", timeout=600)
                print("  [%s] start(attempt %d) rc=%s" % (skill, attempt, rc), flush=True)
                if wait_done(skill, mat):
                    d = collect(skill, mat)
                    print("  [%s] 完成 -> %s" % (skill, json.dumps(d, ensure_ascii=False)[:180]), flush=True)
                    break
                print("  [%s] 未完成，rerun" % skill, flush=True)
                tf_run("-tt", skill, "-p", mat, "-j", "1", "rerun", timeout=120)
            else:
                print("  [%s] 失败" % skill, flush=True)


if __name__ == "__main__":
    main()
