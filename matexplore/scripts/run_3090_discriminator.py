#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""3090 甄别性计算：MatGL 带隙预测(seed vs 替换后候选 对比)。

测试机制 H2/H3：替换是否按预测方向移动 Eg。
输入 POSCAR 来自 /tmp/contrast/(SEED_*.vasp / CAND_*.vasp)，由本地生成。
"""
import csv, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from matexplore.config import load_config

REMOTE_DIR = "/home/wangchaoyue852/matexplore_discriminator"


def ssh(host, cmd, timeout=900):
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, cmd],
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or ""), (r.stderr or "")


def write_remote(host, path, content):
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", host,
                        "mkdir -p $(dirname %s) && cat > %s" % (path, path)],
                       input=content, capture_output=True, text=True, timeout=60)
    return r.returncode


def main():
    cfg = load_config()
    host = cfg.validation.hpc["3090"].ssh_host
    here = os.path.dirname(os.path.abspath(__file__))

    # push matgl script
    write_remote(host, REMOTE_DIR + "/matgl_predict_gap.py",
                 open(os.path.join(here, "matgl_predict_gap.py")).read())
    # push contrast POSCARs
    files = sorted(os.listdir("/tmp/contrast"))
    for fn in files:
        write_remote(host, REMOTE_DIR + "/" + fn, open("/tmp/contrast/" + fn).read())

    paths = " ".join(REMOTE_DIR + "/" + fn for fn in files)
    cmd = "cd %s && ~/matgl_venv/bin/python matgl_predict_gap.py %s" % (REMOTE_DIR, paths)
    print("running MatGL band-gap prediction ...", flush=True)
    rc, out, err = ssh(host, cmd, timeout=900)
    print("exit", rc)
    print("STDOUT:", out[-3000:])
    if err.strip():
        print("STDERR:", err[-3000:])


if __name__ == "__main__":
    main()
