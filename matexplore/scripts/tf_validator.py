#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tf_validator.py —— 经 taskflow 提交做多属性验证链。

计算平台：MACE 链(opt-mace-cpu / phonon-mace-cpu / kl-mace-cpu)走 jzzn(CPU 集群)；
unihamgnn(带隙)走 3090(仅部署在 3090，CPU 计算)。
"""
import argparse, csv, json, os, subprocess, sys

PROJECT_ROOT = "/home/wangchao/software/taskflow/test/tf_test"

# skill -> (summary文件, 摘要步骤目录, hpc, work_dir)
SKILLS = {
    "opt-mace-cpu":    ("energy_summary.json", "step3_formation", "jzzn",
                        "/public/home/wangchao/Fullerene_Network/work"),
    "phonon-mace-cpu": ("phonon_summary.json", "step3_phonon", "jzzn",
                        "/public/home/wangchao/Fullerene_Network/work"),
    "kl-mace-cpu":     ("kappa_summary.json",  "step4_kappa", "jzzn",
                        "/public/home/wangchao/Fullerene_Network/work"),
    "unihamgnn":       ("band_summary.json",   "step3_band", "3090",
                        "/home/wangchaoyue852/taskflow/work"),
}

TF = "/home/wangchao/.local/bin/tf"


def tf_run(*args, timeout=300):
    r = subprocess.run([TF] + list(args), cwd=PROJECT_ROOT,
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or ""), (r.stderr or "")


def setup_material(mat, poscar_text):
    matdir = os.path.join(PROJECT_ROOT, mat)
    os.makedirs(matdir, exist_ok=True)
    with open(os.path.join(matdir, "POSCAR"), "w") as f:
        f.write(poscar_text)
    for skill, (fn, step, hpc, work_dir) in SKILLS.items():
        ps = os.path.join(matdir, skill, "project_setting")
        os.makedirs(ps, exist_ok=True)
        with open(os.path.join(ps, "tf_%s_%s.yaml" % (mat, skill)), "w") as f:
            f.write("task_types:\n  %s:\n    local_root: \"..\"\n"
                    "    work_dir: %s\n" % (skill, work_dir))
        with open(os.path.join(ps, "hpc.yaml"), "w") as f:
            f.write("name: \"%s\"\nssh_host: %s\ntemplate_map: {}\n"
                    % (hpc, hpc))
        with open(os.path.join(ps, "setting.yaml"), "w") as f:
            f.write("auto_advance: false\nbase_dir: \"{matdir}\"\n"
                    "result_dir: \"{matdir}/result\"\nlog_dir: \"{matdir}/log\"\n"
                    "work_dir: %s\n" % work_dir)
    print("setup", mat)


def submit(skill, mat):
    rc, out, err = tf_run("-tt", skill, "-p", mat, "start", timeout=600)
    print("  [%s/%s] rc=%s %s" % (skill, mat, rc, (out + err).strip()[-100:]))
    return rc


def collect(skill, mat):
    fn, step, hpc, work_dir = SKILLS[skill]
    p = os.path.join(PROJECT_ROOT, mat, skill, "result", step, fn)
    if os.path.isfile(p):
        with open(p) as f:
            return json.load(f)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["setup", "submit", "collect", "status"])
    ap.add_argument("--skills", default=",".join(SKILLS.keys()))
    ap.add_argument("--mats", default=None)
    ap.add_argument("--poscar-dir", default="/tmp/contrast2")
    a = ap.parse_args()
    skills = [s for s in a.skills.split(",") if s]

    if a.action == "setup":
        mats = a.mats.split(",") if a.mats else             [fn[:-5] for fn in sorted(os.listdir(a.poscar_dir)) if fn.endswith(".vasp")]
        for mat in mats:
            p = os.path.join(a.poscar_dir, mat + ".vasp")
            if os.path.isfile(p):
                setup_material(mat, open(p).read())
        return

    mats = a.mats.split(",") if a.mats else         [d for d in sorted(os.listdir(PROJECT_ROOT))
         if os.path.isfile(os.path.join(PROJECT_ROOT, d, "POSCAR")) and d != "Si"]

    if a.action == "submit":
        for mat in mats:
            for skill in skills:
                submit(skill, mat)
        return
    if a.action == "status":
        for mat in mats:
            for skill in skills:
                rc, out, err = tf_run("-tt", skill, "-p", mat, "status", timeout=180)
                print("=== %s/%s rc=%s ===" % (skill, mat, rc))
                print((out + err)[:400])
        return
    if a.action == "collect":
        rows = []
        for mat in mats:
            row = {"material": mat}
            for skill in skills:
                d = collect(skill, mat)
                if d is None:
                    row[skill] = ""
                    continue
                row[skill] = "done"
                if skill == "unihamgnn":
                    row["band_gap_eV"] = d.get("band_gap_eV")
                elif skill == "opt-mace-cpu":
                    row["E_per_atom_eV"] = d.get("E_per_atom_eV")
                elif skill == "phonon-mace-cpu":
                    row["phonon_stable"] = d.get("stable")
                elif skill == "kl-mace-cpu":
                    row["kappa"] = d.get("kappa")
            rows.append(row)
        with open(os.path.join(PROJECT_ROOT, "validation_collect.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        for r in rows:
            print(json.dumps(r, ensure_ascii=False))


if __name__ == "__main__":
    main()
