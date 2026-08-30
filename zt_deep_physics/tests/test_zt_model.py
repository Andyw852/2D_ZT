from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from zt_model import (  # noqa: E402
    ElectronicParams,
    LatticeParams,
    StructureParams,
    bipolar_gap_threshold_eV,
    carrier_sheet_density_m2,
    deformation_potential_mobility_m2_Vs,
    electronic_transport,
    lattice_transport,
    spb_seebeck_lorenz,
)


def test_sheet_density_increases_with_eta():
    p = ElectronicParams()
    n = carrier_sheet_density_m2(np.array([-3.0, 0.0, 3.0]), p)
    assert np.all(np.diff(n) > 0)


def test_seebeck_decreases_and_lorenz_increases_toward_degenerate_limit():
    s, lorenz = spb_seebeck_lorenz(np.array([-4.0, 0.0, 6.0]))
    assert np.all(np.diff(s) < 0)
    assert np.all(np.diff(lorenz) > 0)
    assert 2.0e-8 < lorenz[-1] < 2.6e-8


def test_dp_mobility_scalings():
    p = ElectronicParams()
    mu = deformation_potential_mobility_m2_Vs(p)
    assert np.isclose(deformation_potential_mobility_m2_Vs(replace(p, temperature_K=1200)), mu / 2)
    assert np.isclose(
        deformation_potential_mobility_m2_Vs(replace(p, deformation_potential_eV=12)),
        mu / 4,
    )


def test_lattice_velocity_and_overlap_effects():
    lp = LatticeParams()
    base = lattice_transport(lp, StructureParams())
    soft = lattice_transport(lp, StructureParams(stiffness_ratio=0.5))
    overlap = lattice_transport(lp, StructureParams(acoustic_optical_overlap=0.8))
    assert soft["group_velocity_m_s"] < base["group_velocity_m_s"]
    assert soft["kappa_lattice_W_mK"] < base["kappa_lattice_W_mK"]
    assert overlap["lifetime_ps"] < base["lifetime_ps"]
    assert overlap["kappa_lattice_W_mK"] < base["kappa_lattice_W_mK"]


def test_porosity_penalises_both_conductivities_but_zt_remains_finite():
    ep = ElectronicParams()
    lp = LatticeParams()
    dense_s = StructureParams()
    pore_s = StructureParams(porosity=0.25)
    dense_k = lattice_transport(lp, dense_s)["kappa_lattice_W_mK"]
    pore_k = lattice_transport(lp, pore_s)["kappa_lattice_W_mK"]
    dense = electronic_transport(0.0, ep, dense_k, dense_s)
    pore = electronic_transport(0.0, ep, pore_k, pore_s)
    assert pore_k < dense_k
    assert pore["sigma_S_m"] < dense["sigma_S_m"]
    assert np.isfinite(pore["zt"]) and pore["zt"] > 0


def test_gap_threshold_at_600K():
    assert np.isclose(bipolar_gap_threshold_eV(600, 10), 0.51704, rtol=2e-3)

