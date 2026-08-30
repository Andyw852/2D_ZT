#!/usr/bin/env python3
"""快速声子稳定性：不弛豫(输入=已弛豫 CONTCAR)、CPU、2x2x1 超胞。"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from ase import Atoms
from ase.io import read
from mace.calculators import MACECalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

MODEL_DIR = "/home/wangchaoyue852/mace_models"
MODEL = "MACE-matpes-pbe-omat-ft.model"

def ase_to_phonopy(a):
    return PhonopyAtoms(symbols=list(a.get_chemical_symbols()),
                        cell=np.array(a.cell[:], dtype=float),
                        scaled_positions=np.array(a.get_scaled_positions(), dtype=float))

def main(poscar, out):
    atoms = read(poscar, format="vasp")
    calc = MACECalculator(model_paths=[os.path.join(MODEL_DIR, MODEL)],
                          device="cpu", default_dtype="float64")
    unit = ase_to_phonopy(atoms)
    ph = Phonopy(unit, supercell_matrix=np.diag([2, 2, 1]))
    ph.generate_displacements(distance=0.01)
    scs = ph.supercells_with_displacements
    print("n_disp =", len(scs), flush=True)
    forces = []
    for i, sc in enumerate(scs):
        a = Atoms(numbers=sc.numbers, positions=sc.positions, cell=sc.cell, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
        if (i + 1) % 10 == 0:
            print("  force %d/%d" % (i + 1, len(scs)), flush=True)
    ph.forces = np.array(forces, dtype=float)
    ph.produce_force_constants()
    ph.run_mesh([5, 5, 1])
    freqs = ph.get_mesh_dict()["frequencies"]
    fmin = float(freqs.min())
    n_imag = int((freqs < -0.1).sum())
    stable = bool(fmin > -0.3 and n_imag == 0)
    result = {"stable": stable, "min_freq_THz": round(fmin, 4), "n_imag": n_imag,
              "n_atoms": len(atoms), "n_disp": len(scs)}
    print("RESULT:", json.dumps(result, ensure_ascii=False), flush=True)
    with open(out, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
