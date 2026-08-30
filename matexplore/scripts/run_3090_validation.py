#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""推进4：把闭环胜者候选真正提交到 3090 做 MACE 弛豫(T1 稳定性验证)。

复用 taskflow 技能 opt-mace-gpu 的引擎 mace_relax.py + MACE-matpes 模型。
幂等：已标记 T1_mace_relaxed(且 converged)的候选跳过；远程已有 relax_summary.json 则复用。
"""
import argparse, csv, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from matexplore.config import load_config
from matexplore.generation.structure_models import make_poscar

ENGINE_DIR = "/home/wangchao/software/taskflow/skill/_common/mace"
REMOTE_ROOT = "/home/wangchaoyue852/matexplore_validation"
MODEL = "MACE-matpes-pbe-omat-ft.model"
MODEL_DIR = "/home/wangchaoyue852/mace_models"
CUDA_VISIBLE = "1"
PER_CAND_TIMEOUT = 240   # 秒，单个候选的上限


def ssh(host, cmd, timeout=300):
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                        host, cmd], capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or ""), (r.stderr or "")


def write_remote(host, path, content):
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", host,
                        f"mkdir -p $(dirname {path}) && cat > {path}"],
                       input=content, capture_output=True, text=True, timeout=60)
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--formula", default=None)
    ap.add_argument("--force", action="store_true", help="即使已验证也重跑")
    args = ap.parse_args()

    cfg = load_config()
    runs = sorted([d for d in os.listdir(cfg.store.root)
                   if os.path.isdir(os.path.join(cfg.store.root, d))])
    run = runs[-1]
    ledger_path = os.path.join(cfg.store.root, run, "ledger.json")
    ledger = json.load(open(ledger_path))
    winners = ledger["winners"]
    if args.formula:
        winners = [w for w in winners if w["formula"] == args.formula] or winners[:1]
    winners = winners[:args.n]

    host = cfg.validation.hpc["3090"].ssh_host
    done = {k: v for k, v in ledger.get("validated", {}).items()
            if v.get("status") == "T1_mace_relaxed" and v.get("converged")}

    # 推引擎文件(只推一次)
    ssh(host, f"mkdir -p {REMOTE_ROOT}")
    for fn in ("mace_relax.py", "mace_model.py"):
        content = open(os.path.join(ENGINE_DIR, fn)).read()
        assert write_remote(host, f"{REMOTE_ROOT}/{fn}", content) == 0, f"push {fn} failed"

    results = []
    for w in winners:
        formula = w["formula"]
        if not args.force and formula in done:
            print(f"[{formula}] 已 T1_mace_relaxed(converged)，跳过", flush=True)
            continue

        poscar = make_poscar(w["lattice"], w["species"], w["coords"], comment=formula)
        wdir = f"{REMOTE_ROOT}/{formula}"
        write_remote(host, f"{wdir}/POSCAR", poscar + chr(10))
        write_remote(host, f"{wdir}/mace_relax.py", open(os.path.join(ENGINE_DIR, "mace_relax.py")).read())
        write_remote(host, f"{wdir}/mace_model.py", open(os.path.join(ENGINE_DIR, "mace_model.py")).read())

        # 若远程已有 summary，直接复用
        rc_s, summ_txt, _ = ssh(host, f"cat {wdir}/relax_summary.json 2>/dev/null")
        summary = None
        if rc_s == 0 and summ_txt.strip():
            try:
                summary = json.loads(summ_txt)
            except Exception:
                summary = None
        if summary is None:
            cmd = (
                f"cd {wdir} && CUDA_VISIBLE_DEVICES={CUDA_VISIBLE} "
                f"bash -lc 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate mace-gpu && "
                f"python mace_relax.py --model {MODEL} --model-dir {MODEL_DIR} --dim 2d --device cuda "
                f"--cell-policy none --relax true --relax-cell true --fmax 1e-3 --steps 600 "
                f"--fix-symmetry true --residual-tol 2e-3 --stress-tol 0.5'"
            )
            print(f"[{formula}] 提交 MACE 弛豫 ...", flush=True)
            try:
                rc, out, err = ssh(host, cmd, timeout=PER_CAND_TIMEOUT)
            except subprocess.TimeoutExpired:
                rc, out, err = -1, "", f"TIMEOUT>{PER_CAND_TIMEOUT}s"
            rc_s, summ_txt, _ = ssh(host, f"cat {wdir}/relax_summary.json 2>/dev/null")
            summary = json.loads(summ_txt) if (rc_s == 0 and summ_txt.strip()) else None

        row = {"formula": formula, "seed": w.get("seed")}
        if summary:
            row.update({
                "converged": summary.get("converged"),
                "E_eV": round(summary.get("energy_eV", 0), 4),
                "E_per_atom_eV": round(summary.get("energy_per_atom_eV", 0), 4),
                "max_force_eV_A": round(summary.get("max_force_eV_per_A", 0), 6),
                "max_stress_GPa": round(summary.get("max_stress_gated_GPa", 0), 6),
                "vol_change_pct": round(summary.get("volume_change_pct", 0), 3),
                "spg_in": summary.get("spacegroup_in"),
                "spg_out": summary.get("spacegroup_out"),
                "steps": summary.get("opt_steps"),
                "wall_s": summary.get("wall_time_s"),
            })
        else:
            row.update({"converged": False, "note": "no summary", "tail": (out + err)[-300:]})
        results.append(row)
        print("   ->", json.dumps({k: row.get(k) for k in
              ("converged", "E_per_atom_eV", "max_force_eV_A", "max_stress_GPa", "spg_out")},
              ensure_ascii=False, default=str), flush=True)

    if results:
        out_csv = os.path.join(cfg.store.root, run, "validation_3090.csv")
        with open(out_csv, "w", newline="") as f:
            wcsv = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            wcsv.writeheader()
            wcsv.writerows(results)
        for row in results:
            ledger.setdefault("validated", {})[row["formula"]] = {
                "status": "T1_mace_relaxed",
                "converged": row.get("converged"),
                "E_per_atom_eV": row.get("E_per_atom_eV"),
                "max_force_eV_A": row.get("max_force_eV_A"),
                "spg_out": row.get("spg_out"),
            }
        json.dump(ledger, open(ledger_path, "w"), ensure_ascii=False, indent=2, default=str)
        print("\n写回:", out_csv)

    # 汇总打印
    print("\n=== 3090 MACE 弛豫结果 ===")
    for row in results:
        print(f"  {row['formula']:14s} converged={row.get('converged')}  "
              f"E/atom={row.get('E_per_atom_eV')} eV  "
              f"maxF={row.get('max_force_eV_A')}  stress={row.get('max_stress_GPa')} GPa")


if __name__ == "__main__":
    main()

