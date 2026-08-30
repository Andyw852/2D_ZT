#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""推进 4 属性验证链直到全部完成(监控循环，后台跑)。"""
import os, sys, time, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tf_validator as tv

PROJECT_ROOT = tv.PROJECT_ROOT
SKILLS = tv.SKILLS
tf_run = tv.tf_run
collect = tv.collect

MATS = [d for d in sorted(os.listdir(PROJECT_ROOT))
        if os.path.isfile(os.path.join(PROJECT_ROOT, d, "POSCAR"))
        and d not in ("Si",)]


def main():
    skills = list(SKILLS.keys())
    total = len(MATS) * len(skills)
    print("监控循环开始: %d 材料 x %d 技能 = %d 项" % (len(MATS), len(skills), total), flush=True)
    for it in range(1, 400):
        for skill in skills:
            rc, out, err = tf_run("-tt", skill, "start", timeout=400)
            n_sub = (out + err).count("已提交")
            if n_sub:
                print("  [%s] 本轮提交 %d 个作业" % (skill, n_sub), flush=True)
        done = 0
        for mat in MATS:
            for skill in skills:
                if collect(skill, mat) is not None:
                    done += 1
        print("[iter %d] %d/%d 完成" % (it, done, total), flush=True)
        if done >= total:
            print("全部完成", flush=True)
            break
        time.sleep(150)


if __name__ == "__main__":
    main()
