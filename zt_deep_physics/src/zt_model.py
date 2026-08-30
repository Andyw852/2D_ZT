"""Transparent reduced-order models for 2D thermoelectric design.

The electronic model is a dimensionally consistent single-parabolic-band
model.  The lattice model is deliberately a reduced modal-BTE scenario model;
it is useful for sensitivity and design-window studies, not material-specific
prediction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Any

import numpy as np
from scipy.integrate import quad
from scipy.special import expit


KB = 1.380649e-23
E_CHARGE = 1.602176634e-19
HBAR = 1.054571817e-34
M_E = 9.1093837015e-31


@dataclass(frozen=True)
class ElectronicParams:
    temperature_K: float = 600.0
    dos_mass_me: float = 0.8
    conductivity_mass_me: float = 0.35
    valley_degeneracy: int = 2
    elastic_modulus_2d_N_per_m: float = 50.0
    deformation_potential_eV: float = 6.0
    effective_thickness_nm: float = 1.0
    scattering_exponent: float = 0.0
    mobility_retention: float = 0.85


@dataclass(frozen=True)
class LatticeParams:
    heat_capacity_J_per_m3K: float = 1.6e6
    sound_velocity_m_per_s: float = 2500.0
    intrinsic_lifetime_ps_at_600K: float = 0.50
    temperature_K: float = 600.0


@dataclass(frozen=True)
class StructureParams:
    porosity: float = 0.0
    pore_spacing_nm: float = 20.0
    wrinkle_slope: float = 0.0
    wrinkle_wavelength_nm: float = 20.0
    stiffness_ratio: float = 1.0
    mass_ratio: float = 1.0
    acoustic_optical_overlap: float = 0.15


def _validate_positive(name: str, value: float) -> None:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and > 0, got {value!r}")


@lru_cache(maxsize=16384)
def _fermi_scalar(order: float, eta: float) -> float:
    """Unnormalised complete Fermi-Dirac integral F_order(eta)."""
    if order <= -1:
        raise ValueError("Fermi-integral order must be > -1")
    if abs(order) < 1e-14:
        return float(np.logaddexp(0.0, eta))
    value, _ = quad(
        lambda x: x**order * expit(eta - x),
        0.0,
        np.inf,
        epsabs=2e-10,
        epsrel=2e-10,
        limit=250,
    )
    return float(value)


def fermi_integral(order: float, eta: Any) -> np.ndarray | float:
    arr = np.asarray(eta, dtype=float)
    flat = np.array([_fermi_scalar(float(order), float(v)) for v in arr.ravel()])
    out = flat.reshape(arr.shape)
    return float(out) if out.ndim == 0 else out


def spb_seebeck_lorenz(
    eta: Any,
    dimension: int = 2,
    scattering_exponent: float = 0.0,
) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Return |S| [V/K] and Lorenz number [W ohm K^-2].

    For a parabolic band the transport distribution is proportional to E^a,
    where a = r + d/2.  This implementation uses unnormalised Fermi integrals.
    """
    if dimension not in (2, 3):
        raise ValueError("dimension must be 2 or 3")
    a = scattering_exponent + dimension / 2.0
    if a <= 0:
        raise ValueError("r + d/2 must be > 0")
    eta_arr = np.asarray(eta, dtype=float)
    f0 = np.asarray(fermi_integral(a - 1.0, eta_arr))
    f1 = np.asarray(fermi_integral(a, eta_arr))
    f2 = np.asarray(fermi_integral(a + 1.0, eta_arr))
    mean_reduced_energy = (a + 1.0) * f1 / (a * f0)
    seebeck = (KB / E_CHARGE) * (mean_reduced_energy - eta_arr)
    lorenz = (KB / E_CHARGE) ** 2 * (
        (a + 2.0) * f2 / (a * f0) - mean_reduced_energy**2
    )
    if eta_arr.ndim == 0:
        return float(seebeck), float(lorenz)
    return seebeck, lorenz


def carrier_sheet_density_m2(eta: Any, p: ElectronicParams) -> np.ndarray | float:
    """2D sheet density including spin degeneracy and explicit valley count."""
    _validate_positive("temperature_K", p.temperature_K)
    _validate_positive("dos_mass_me", p.dos_mass_me)
    if p.valley_degeneracy < 1:
        raise ValueError("valley_degeneracy must be >= 1")
    prefactor = (
        p.valley_degeneracy
        * p.dos_mass_me
        * M_E
        * KB
        * p.temperature_K
        / (np.pi * HBAR**2)
    )
    return prefactor * fermi_integral(0.0, eta)


def deformation_potential_mobility_m2_Vs(p: ElectronicParams) -> float:
    """Acoustic-deformation-potential mobility for a 2D parabolic band."""
    for name in (
        "temperature_K",
        "dos_mass_me",
        "conductivity_mass_me",
        "elastic_modulus_2d_N_per_m",
        "deformation_potential_eV",
    ):
        _validate_positive(name, float(getattr(p, name)))
    if not (0 < p.mobility_retention <= 1):
        raise ValueError("mobility_retention must be in (0, 1]")
    md = p.dos_mass_me * M_E
    mc = p.conductivity_mass_me * M_E
    e1 = p.deformation_potential_eV * E_CHARGE
    mobility = (
        E_CHARGE
        * HBAR**3
        * p.elastic_modulus_2d_N_per_m
        / (KB * p.temperature_K * mc * md * e1**2)
    )
    return float(mobility * p.mobility_retention)


def structure_electrical_factor(s: StructureParams) -> float:
    """Scenario-level retention of connected electronic transport paths."""
    if not (0 <= s.porosity < 0.75):
        raise ValueError("porosity must be in [0, 0.75)")
    _validate_positive("stiffness_ratio", s.stiffness_ratio)
    tortuosity = 1.0 + 0.25 * s.wrinkle_slope**2
    pore_connectivity = (1.0 - s.porosity) ** 1.5
    curvature_retention = 1.0 / (1.0 + 0.15 * s.wrinkle_slope**2)
    # In the deformation-potential approximation mu is linear in C_2D.  The
    # stiffness ratio therefore acts on the reference mobility as well as on
    # the phonon velocity.  This is the central soft-lattice trade-off.
    return float(s.stiffness_ratio * pore_connectivity * curvature_retention / tortuosity)


def lattice_transport(p: LatticeParams, s: StructureParams) -> dict[str, float]:
    """Reduced modal-BTE lattice conductivity with transparent modifiers."""
    for name in (
        "heat_capacity_J_per_m3K",
        "sound_velocity_m_per_s",
        "intrinsic_lifetime_ps_at_600K",
        "temperature_K",
    ):
        _validate_positive(name, float(getattr(p, name)))
    for name in ("pore_spacing_nm", "wrinkle_wavelength_nm", "stiffness_ratio", "mass_ratio"):
        _validate_positive(name, float(getattr(s, name)))
    if not (0 <= s.porosity < 0.75):
        raise ValueError("porosity must be in [0, 0.75)")
    if not (0 <= s.acoustic_optical_overlap <= 1):
        raise ValueError("acoustic_optical_overlap must be in [0, 1]")

    velocity = p.sound_velocity_m_per_s * np.sqrt(s.stiffness_ratio / s.mass_ratio)
    tortuosity = 1.0 + 0.25 * s.wrinkle_slope**2
    projected_velocity = velocity / tortuosity

    tau_intrinsic = (
        p.intrinsic_lifetime_ps_at_600K
        * 1e-12
        * 600.0
        / p.temperature_K
        / (1.0 + 3.0 * s.acoustic_optical_overlap)
    )
    pore_rate = velocity * s.porosity / (s.pore_spacing_nm * 1e-9)
    wrinkle_rate = velocity * s.wrinkle_slope**2 / (s.wrinkle_wavelength_nm * 1e-9)
    tau_effective = 1.0 / (1.0 / tau_intrinsic + pore_rate + wrinkle_rate)

    effective_medium = (1.0 - s.porosity) / (1.0 + s.porosity)
    kappa_lattice = (
        0.5
        * p.heat_capacity_J_per_m3K
        * projected_velocity**2
        * tau_effective
        * effective_medium
    )
    return {
        "kappa_lattice_W_mK": float(kappa_lattice),
        "group_velocity_m_s": float(projected_velocity),
        "lifetime_ps": float(tau_effective * 1e12),
        "mean_free_path_nm": float(projected_velocity * tau_effective * 1e9),
        "tortuosity": float(tortuosity),
        "effective_medium_factor": float(effective_medium),
    }


def electronic_transport(
    eta: Any,
    p: ElectronicParams,
    kappa_lattice_W_mK: float,
    structure: StructureParams | None = None,
) -> dict[str, np.ndarray | float]:
    """Return electronic transport and ZT for scalar or array eta."""
    _validate_positive("effective_thickness_nm", p.effective_thickness_nm)
    if kappa_lattice_W_mK < 0:
        raise ValueError("kappa_lattice_W_mK must be >= 0")
    sheet_density = np.asarray(carrier_sheet_density_m2(eta, p))
    thickness_m = p.effective_thickness_nm * 1e-9
    volume_density = sheet_density / thickness_m
    mobility = deformation_potential_mobility_m2_Vs(p)
    electrical_factor = structure_electrical_factor(structure) if structure else 1.0
    sigma = volume_density * E_CHARGE * mobility * electrical_factor
    seebeck, lorenz = spb_seebeck_lorenz(
        eta, dimension=2, scattering_exponent=p.scattering_exponent
    )
    seebeck = np.asarray(seebeck)
    lorenz = np.asarray(lorenz)
    power_factor = seebeck**2 * sigma
    kappa_e = lorenz * sigma * p.temperature_K
    zt = power_factor * p.temperature_K / (kappa_e + kappa_lattice_W_mK)
    result = {
        "eta": np.asarray(eta, dtype=float),
        "n_sheet_cm2": sheet_density / 1e4,
        "n_equiv_cm3": volume_density / 1e6,
        "seebeck_uV_K": seebeck * 1e6,
        "mobility_cm2_Vs": mobility * 1e4 * electrical_factor,
        "sigma_S_m": sigma,
        "lorenz_W_ohm_K2": lorenz,
        "power_factor_mW_mK2": power_factor * 1e3,
        "kappa_e_W_mK": kappa_e,
        "kappa_lattice_W_mK": np.zeros_like(sigma) + kappa_lattice_W_mK,
        "zt": zt,
    }
    if np.asarray(eta).ndim == 0:
        return {key: float(np.asarray(value)) for key, value in result.items()}
    return result


def optimise_eta(
    eta_grid: np.ndarray,
    p: ElectronicParams,
    lattice: LatticeParams,
    structure: StructureParams,
) -> dict[str, float]:
    lp = replace(lattice, temperature_K=p.temperature_K)
    k_l = lattice_transport(lp, structure)["kappa_lattice_W_mK"]
    values = electronic_transport(eta_grid, p, k_l, structure)
    idx = int(np.nanargmax(values["zt"]))
    result = {}
    for key, value in values.items():
        arr = np.asarray(value)
        result[key] = float(arr) if arr.ndim == 0 else float(arr[idx])
    return result


def bipolar_gap_threshold_eV(temperature_K: float, multiple_kBT: float = 10.0) -> float:
    """Rule-of-thumb minimum gap, not an upper or universal optimum."""
    _validate_positive("temperature_K", temperature_K)
    _validate_positive("multiple_kBT", multiple_kBT)
    return float(multiple_kBT * KB * temperature_K / E_CHARGE)
