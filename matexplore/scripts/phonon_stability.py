#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""声子动力学稳定性(MACE + phonopy)：MACE 弛豫 -> 超胞位移 -> fc2 -> 虚频闸。

在 3090 mace-gpu 环境跑。用法: python phonon_stability.py POSCAR out.json
输出: {"stable": bool, "min_freq_THz": float, "n_imag": int, "relax_converged": bool}
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from ase import Atoms
from ase.io import read
from ase.optimize import FIRE
from mace.calculators import MACECalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

MODEL_DIR = "/home/wangchaoyue852/mace_models"
MODEL = "MACE-matpes-pbe-omat-ft.model"


def ase_to_phonopy(a):
    return PhonopyAtoms(symbols=list(a.get_chemical_symbols()),
                        cell=np.array(a.cell[:], dtype=float),
                        scaled_positions=np.array(a.get_scaled_positions(), dtype=float))


def relax_2d(atoms, calc):
    """2D 弛豫：放开面内轴长+面内剪切，锁死真空方向(Voigt [xx,yy,zz,yz,xz,xy])。"""
    from ase.filters import FrechetCellFilter
    mask = [1, 1, 0, 0, 0, 1]  # 锁 zz,yz,xz
    atoms.calc = calc
    cf = FrechetCellFilter(atoms, mask=mask)
    opt = FIRE(cf, logfile=None)
    try:
        opt.run(fmax=1e-3, steps=400)
        return atoms, bool(np.max(np.linalg.norm(cf.get_forces(), axis=1)) < 1e-3)
    except Exception as e:
        return atoms, False


def main(poscar, out):
    atoms = read(poscar, format="vasp")
    calc = MACECalculator(model_paths=[os.path.join(MODEL_DIR, MODEL)],
                          device="cuda", default_dtype="float64")
    atoms, relaxed = relax_2d(atoms, calc)
    atoms.calc = calc

    unit = ase_to_phonopy(atoms)
    ph = Phonopy(unit, supercell_matrix=np.diag([3, 3, 1]))
    ph.generate_displacements(distance=0.01)
    scs = ph.supercells_with_displacements
    forces = []
    for sc in scs:
        a = Atoms(numbers=sc.numbers, positions=sc.positions, cell=sc.cell, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
    ph.forces = np.array(forces, dtype=float)
    ph.produce_force_constants()
    ph.run_mesh([11, 11, 1])
    freqs = ph.get_mesh_dict()["frequencies"]
    fmin = float(freqs.min())
    n_imag = int((freqs < -0.1).sum())
    stable = bool(fmin > -0.3 and n_imag == 0)
    result = {"stable": stable, "min_freq_THz": round(fmin, 4),
              "n_imag": n_imag, "relax_converged": relaxed,
              "n_atoms": len(atoms), "n_disp": len(scs)}
    print(json.dumps(result, ensure_ascii=False))
    with open(out, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
