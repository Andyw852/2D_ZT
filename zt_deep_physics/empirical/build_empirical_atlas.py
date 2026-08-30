"""Build empirical, maximum-coverage thermoelectric relationship figures.

Every panel uses pairwise-complete observations from the data source relevant
to that relation.  Real-ZT panels use one peak-ZT record per experimental
sample and interpolate other properties only within the same sample.
"""

from __future__ import annotations

import ast
import json
import os
import re
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/zt_empirical_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
REPO = PROJECT.parent
STAR = REPO / "jarvis_2d_te_atlas" / "data" / "raw" / "external" / "starrydata2"
OUT = HERE / "outputs"
FIG = HERE / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


PROP_MAP = {
    "ZT": "zt",
    "Seebeck coefficient": "seebeck",
    "Electrical conductivity": "sigma",
    "Electrical resistivity": "resistivity",
    "Power factor": "power_factor",
    "Lattice thermal conductivity": "kappa_lattice",
    "Thermal conductivity": "kappa_total",
    "total thermal conductivity": "kappa_total",
    "Electronic thermal conductivity": "kappa_electronic",
    "electron thermal conductivity": "kappa_electronic",
    "Carrier concentration": "carrier_concentration",
    "Carrier mobility": "mobility",
    "Hall mobility": "mobility",
}

VALID_RANGES = {
    "zt": (1e-8, 5.0),
    "seebeck": (-5e-3, 5e-3),
    "sigma": (1e-3, 1e9),
    "resistivity": (1e-10, 1e3),
    "power_factor": (0.0, 0.2),
    "kappa_lattice": (1e-3, 200.0),
    "kappa_total": (1e-3, 300.0),
    "kappa_electronic": (0.0, 200.0),
    "carrier_concentration": (1e20, 1e30),
    "mobility": (1e-8, 10.0),
}

DISPLAY = {
    "zt": "ZT",
    "seebeck": "Seebeck",
    "sigma": "electrical conductivity",
    "power_factor": "power factor",
    "kappa_lattice": "lattice thermal conductivity",
    "kappa_total": "total thermal conductivity",
    "kappa_electronic": "electronic thermal conductivity",
    "carrier_concentration": "carrier concentration",
    "mobility": "mobility",
}


def parse_array(value) -> np.ndarray:
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
        return np.atleast_1d(np.asarray(parsed, dtype=float))
    except Exception:
        return np.array([], dtype=float)


def clean_curve(prop: str, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x) & np.isfinite(y) & (x >= 100.0) & (x <= 1500.0)
    lo, hi = VALID_RANGES[prop]
    valid &= (y >= lo) & (y <= hi)
    x, y = x[valid], y[valid]
    if len(x) == 0:
        return x, y
    order = np.argsort(x)
    x, y = x[order], y[order]
    unique = np.unique(x)
    y_median = np.array([np.median(y[x == value]) for value in unique])
    return unique, y_median


def load_curve_store():
    columns = ["sample_id", "composition", "SID", "DOI", "prop_x", "prop_y", "unit_y", "x", "y"]
    data = pd.read_csv(STAR / "ThermoelectricMaterials_curves.csv", usecols=columns)
    data = data[data.prop_x.eq("Temperature") & data.prop_y.isin(PROP_MAP)].copy()
    curves: dict[tuple[int, str], list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
    metadata = {}
    for row in data.itertuples(index=False):
        x, y = parse_array(row.x), parse_array(row.y)
        if len(x) == 0 or len(x) != len(y):
            continue
        prop = PROP_MAP[row.prop_y]
        x, y = clean_curve(prop, x, y)
        if len(x) == 0:
            continue
        sid = int(row.sample_id)
        curves[(sid, prop)].append((x, y))
        metadata[sid] = {
            "composition": row.composition,
            "paper_SID": row.SID,
            "DOI": row.DOI,
        }
    return curves, metadata


def values_at_temperature(curves, sample_id: int, prop: str, temperature: float, max_gap=100.0):
    values = []
    for x, y in curves.get((sample_id, prop), []):
        if x.min() <= temperature <= x.max() and np.min(np.abs(x - temperature)) <= max_gap:
            values.append(float(np.interp(temperature, x, y)))
    return values


def median_at_temperature(curves, sample_id: int, prop: str, temperature: float):
    values = values_at_temperature(curves, sample_id, prop, temperature)
    if prop == "sigma" and not values:
        resistivity = values_at_temperature(curves, sample_id, "resistivity", temperature)
        values = [1.0 / value for value in resistivity if value > 0]
    return float(np.median(values)) if values else np.nan


def add_derived_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["abs_seebeck_uV_K"] = frame.seebeck.abs() * 1e6
    frame["carrier_concentration_cm3"] = frame.carrier_concentration / 1e6
    frame["mobility_cm2_Vs"] = frame.mobility * 1e4
    frame["power_factor_mW_mK2"] = frame.power_factor * 1e3
    derived_pf = frame.seebeck**2 * frame.sigma
    frame["power_factor_used_W_mK2"] = frame.power_factor.where(frame.power_factor.notna(), derived_pf)
    frame["power_factor_used_mW_mK2"] = frame.power_factor_used_W_mK2 * 1e3
    return frame


def build_peak_table(curves, metadata):
    properties = [
        "seebeck",
        "sigma",
        "power_factor",
        "kappa_lattice",
        "kappa_total",
        "kappa_electronic",
        "carrier_concentration",
        "mobility",
    ]
    rows = []
    for (sample_id, prop), sample_curves in curves.items():
        if prop != "zt":
            continue
        zt, temperature = max(
            (float(y.max()), float(x[np.argmax(y)])) for x, y in sample_curves
        )
        row = {
            "sample_id": sample_id,
            "zt_peak": zt,
            "temperature_peak_K": temperature,
            **metadata.get(sample_id, {}),
        }
        for target in properties:
            row[target] = median_at_temperature(curves, sample_id, target, temperature)
        rows.append(row)
    frame = add_derived_columns(pd.DataFrame(rows))
    frame.to_csv(OUT / "experimental_peak_ZT_pairs.csv", index=False)
    return frame


def build_common_temperature_table(curves, metadata, temperatures=(300.0, 600.0, 900.0)):
    properties = [
        "zt",
        "seebeck",
        "sigma",
        "power_factor",
        "kappa_lattice",
        "kappa_total",
        "kappa_electronic",
        "carrier_concentration",
        "mobility",
    ]
    sample_ids = sorted({sample_id for sample_id, _ in curves})
    rows = []
    for temperature in temperatures:
        for sample_id in sample_ids:
            row = {"sample_id": sample_id, "temperature_K": temperature, **metadata.get(sample_id, {})}
            any_value = False
            for prop in properties:
                value = median_at_temperature(curves, sample_id, prop, temperature)
                row[prop] = value
                any_value |= np.isfinite(value)
            if any_value:
                rows.append(row)
    frame = add_derived_columns(pd.DataFrame(rows))
    frame.to_csv(OUT / "common_temperature_pairs.csv", index=False)
    return frame


def build_kappa_temperature_band(curves):
    sample_ids = sorted({sample_id for sample_id, prop in curves if prop == "kappa_lattice"})
    rows = []
    for temperature in np.arange(250.0, 1100.1, 25.0):
        for sample_id in sample_ids:
            value = median_at_temperature(curves, sample_id, "kappa_lattice", temperature)
            if np.isfinite(value):
                rows.append({"sample_id": sample_id, "temperature_K": temperature, "kappa_lattice": value})
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "kappa_lattice_temperature_grid.csv", index=False)
    return frame


def normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def info_field(obj, target: str):
    for key, value in obj.items():
        if normalized_key(key) == target and isinstance(value, dict):
            return str(value.get("category", "")).strip(), str(value.get("comment", "")).strip()
    return "", ""


def parse_relative_density(category: str, comment: str) -> float:
    text = f"{category} {comment}".lower()
    if not text.strip() or "unknown" in text:
        return np.nan
    values = [float(value) for value in re.findall(r"(?<![a-z])([0-9]+(?:\.[0-9]+)?)\s*%?", text)]
    values = [value for value in values if 20 <= value <= 100]
    return float(np.median(values)) if values else np.nan


def parse_grain_size_um(category: str, comment: str) -> float:
    # A free-text comment commonly contains the measured range, whereas the
    # category is only a coarse bin.  Do not average the two representations.
    source = comment if comment.strip() else category
    text = source.lower().replace("μm", "um").replace("µm", "um")
    text = text.replace("micrometers", "um").replace("micrometer", "um")
    text = text.replace("nanometers", "nm").replace("nanometer", "nm")
    text = re.sub(r"\*\s*10\^?\(?-6\)?\s*m", " um", text)
    matches = re.findall(
        r"([0-9]+(?:\.[0-9]+)?)\s*(?:to|[-~])?\s*([0-9]+(?:\.[0-9]+)?)?\s*(nm|um|mm)",
        text,
    )
    values = []
    scale = {"nm": 1e-3, "um": 1.0, "mm": 1e3}
    for first, second, unit in matches:
        pair = [float(first)] + ([float(second)] if second else [])
        values.append(float(np.mean(pair)) * scale[unit])
    values = [value for value in values if 0.003 <= value <= 500]
    return float(np.median(values)) if values else np.nan


def canonical_form(value: str) -> str:
    key = normalized_key(value)
    mapping = {
        "bulk": "Bulk",
        "singlecrystal": "Single crystal",
        "polycrystal": "Polycrystal",
        "film": "Film",
        "epitaxialfilm": "Film",
        "multilayerfilm": "Film",
        "pellets": "Pellet/compact",
        "compact": "Pellet/compact",
        "powder": "Powder",
        "wire": "Wire/fiber",
        "fibers": "Wire/fiber",
        "ribbon": "Ribbon/sheet",
        "sheet": "Ribbon/sheet",
        "aerogel": "Porous/aerogel",
        "foamedbulk": "Porous/aerogel",
        "orientedbulk": "Oriented bulk",
    }
    return mapping.get(key, value.strip() if value and value.strip() not in {"-", "Unknown"} else np.nan)


def load_sample_metadata():
    samples = pd.read_csv(
        STAR / "ThermoelectricMaterials_samples.csv",
        usecols=["sample_id", "sample_info"],
    )
    rows = []
    for row in samples.itertuples(index=False):
        try:
            obj = json.loads(row.sample_info) if isinstance(row.sample_info, str) else {}
        except Exception:
            obj = {}
        density_cat, density_comment = info_field(obj, "relativedensity")
        grain_cat, grain_comment = info_field(obj, "grainsize")
        form_cat, _ = info_field(obj, "form")
        family_cat, _ = info_field(obj, "materialfamily")
        relative_density = parse_relative_density(density_cat, density_comment)
        grain_um = parse_grain_size_um(grain_cat, grain_comment)
        rows.append(
            {
                "sample_id": int(row.sample_id),
                "relative_density_pct": relative_density,
                "porosity_fraction": 1.0 - relative_density / 100.0 if np.isfinite(relative_density) else np.nan,
                "grain_size_um": grain_um,
                "form": canonical_form(form_cat),
                "material_family": family_cat if family_cat else np.nan,
            }
        )
    frame = pd.DataFrame(rows).drop_duplicates("sample_id", keep="last")
    frame.to_csv(OUT / "sample_structure_metadata.csv", index=False)
    return frame


def load_elastic_kappa_crossmatch():
    path = REPO / "jarvis_2d_te_atlas" / "features" / "kl_verify" / "kl_views.parquet"
    frame = pd.read_parquet(path)
    for column in ["B_kv", "G_gv", "density", "kL_300"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame[(frame.B_kv > 0) & (frame.G_gv > 0) & (frame.density > 0) & (frame.kL_300 > 0)].copy()
    frame["bulk_sound_speed_proxy_m_s"] = np.sqrt(frame.B_kv * 1e9 / (frame.density * 1000.0))
    frame["shear_sound_speed_proxy_m_s"] = np.sqrt(frame.G_gv * 1e9 / (frame.density * 1000.0))
    frame.to_csv(OUT / "elastic_kappa_crossmatch.csv", index=False)
    return frame


def as_list(value):
    return ast.literal_eval(value) if isinstance(value, str) else value


def build_jarvis2d_geometry():
    structures = pd.read_parquet(
        REPO / "jarvis_2d_te_atlas" / "data" / "processed" / "standardized_2d_structures.parquet"
    )
    rows = []
    for row in structures.itertuples(index=False):
        lattice = np.asarray(as_list(row.lattice), dtype=float)
        positions = np.asarray(as_list(row.positions), dtype=float)
        area_vector = np.cross(lattice[0], lattice[1])
        area = np.linalg.norm(area_vector)
        normal = area_vector / area
        projected = positions @ normal
        thickness = float(projected.max() - projected.min()) if len(projected) else np.nan
        rows.append(
            {
                "jid": row.jid,
                "formula": row.formula,
                "nsites": row.nsites,
                "area_A2": area,
                "area_per_atom_A2": area / row.nsites,
                "out_of_plane_span_A": thickness,
                "out_of_plane_ratio": thickness / np.sqrt(area),
            }
        )
    geometry = pd.DataFrame(rows)
    transport = pd.read_csv(
        REPO / "jarvis_2d_te_atlas" / "data" / "processed" / "ZT_e_all.csv"
    )
    frame = transport.merge(geometry, on=["jid", "formula"], how="inner")
    frame["abs_seebeck_uV_K"] = frame.S_median.abs()
    frame["sigma_S_m"] = 10 ** frame.log_sigma_dom_geo
    frame.to_csv(OUT / "jarvis2d_geometry_transport.csv", index=False)
    return frame


def finite_pairs(frame, x, y, positive_x=False, positive_y=False):
    subset = frame.dropna(subset=[x, y]).copy()
    subset = subset[np.isfinite(subset[x]) & np.isfinite(subset[y])]
    if positive_x:
        subset = subset[subset[x] > 0]
    if positive_y:
        subset = subset[subset[y] > 0]
    return subset


def relation_stats(frame, x, y):
    subset = finite_pairs(frame, x, y)
    if len(subset) < 3:
        return len(subset), np.nan
    return len(subset), float(spearmanr(subset[x], subset[y]).statistic)


def add_stats(ax, frame, x, y, unique_sample_col=None):
    n, rho = relation_stats(frame, x, y)
    sample_text = ""
    if unique_sample_col and unique_sample_col in frame:
        sample_text = f", samples={frame.dropna(subset=[x,y])[unique_sample_col].nunique()}"
    ax.text(
        0.03,
        0.97,
        f"N={n}{sample_text}\nSpearman rho={rho:+.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox=dict(fc="white", ec="0.7", alpha=0.85),
    )


def add_binned_summary(ax, x, y, log_x=False, bins=12, color="white"):
    x, y = np.asarray(x, float), np.asarray(y, float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0 if log_x else True)
    x, y = x[valid], y[valid]
    if len(x) < 30:
        return
    lo, hi = np.quantile(x, [0.01, 0.99])
    edges = np.geomspace(lo, hi, bins + 1) if log_x else np.linspace(lo, hi, bins + 1)
    centers, medians, lower, upper = [], [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        values = y[(x >= left) & (x < right)]
        if len(values) < 10:
            continue
        centers.append(np.sqrt(left * right) if log_x else 0.5 * (left + right))
        medians.append(np.median(values))
        lower.append(np.quantile(values, 0.25))
        upper.append(np.quantile(values, 0.75))
    if centers:
        ax.plot(centers, medians, color=color, lw=2.0, marker="o", ms=3, label="_nolegend_")
        ax.fill_between(centers, lower, upper, color=color, alpha=0.18, label="_nolegend_")


def density_panel(ax, frame, x, y, xlabel, ylabel, log_x=False, log_y=False, title=""):
    subset = finite_pairs(frame, x, y, positive_x=log_x, positive_y=log_y)
    if log_x:
        x_plot = np.log10(subset[x])
        xlabel = "log10 " + xlabel
    else:
        x_plot = subset[x]
    if log_y:
        y_plot = np.log10(subset[y])
        ylabel = "log10 " + ylabel
    else:
        y_plot = subset[y]
    hb = ax.hexbin(x_plot, y_plot, gridsize=45, mincnt=1, bins="log", cmap="viridis")
    add_binned_summary(ax, x_plot, y_plot, log_x=False, bins=12, color="white")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    add_stats(ax, subset.assign(_x=x_plot, _y=y_plot), "_x", "_y", "sample_id")
    return hb, len(subset)


def make_coverage_figure(curves, peak, common, metadata, elastic, jarvis2d):
    direct_counts = {}
    for prop in DISPLAY:
        direct_counts[prop] = len({sample_id for sample_id, p in curves if p == prop})
    peak_pairs = {
        "ZT + n": peak[["zt_peak", "carrier_concentration"]].notna().all(axis=1).sum(),
        "ZT + S": peak[["zt_peak", "seebeck"]].notna().all(axis=1).sum(),
        "ZT + sigma": peak[["zt_peak", "sigma"]].notna().all(axis=1).sum(),
        "ZT + PF": peak[["zt_peak", "power_factor_used_W_mK2"]].notna().all(axis=1).sum(),
        "ZT + kL": peak[["zt_peak", "kappa_lattice"]].notna().all(axis=1).sum(),
        "ZT + kTotal": peak[["zt_peak", "kappa_total"]].notna().all(axis=1).sum(),
    }
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.8))
    labels = [DISPLAY[key] for key in direct_counts]
    values = list(direct_counts.values())
    order = np.argsort(values)
    axs[0, 0].barh(np.array(labels)[order], np.array(values)[order], color="#4b8bc4")
    for i, value in enumerate(np.array(values)[order]):
        axs[0, 0].text(value, i, f" {value:,}", va="center", fontsize=8)
    axs[0, 0].set_xlabel("unique experimental samples")
    axs[0, 0].set_title("(a) Starrydata2 direct property coverage")

    pair_labels, pair_values = list(peak_pairs), list(peak_pairs.values())
    order = np.argsort(pair_values)
    axs[0, 1].barh(np.array(pair_labels)[order], np.array(pair_values)[order], color="#d47a42")
    for i, value in enumerate(np.array(pair_values)[order]):
        axs[0, 1].text(value, i, f" {value:,}", va="center", fontsize=8)
    axs[0, 1].set_xlabel("same-sample pairs at peak-ZT temperature")
    axs[0, 1].set_title("(b) Real-ZT panel sample sizes")

    temps = sorted(common.temperature_K.unique())
    properties = ["zt", "seebeck", "sigma", "kappa_lattice", "carrier_concentration"]
    matrix = np.array(
        [[common.loc[common.temperature_K.eq(temp), prop].notna().sum() for prop in properties] for temp in temps]
    )
    im = axs[1, 0].imshow(np.log10(matrix + 1), cmap="Blues", aspect="auto")
    axs[1, 0].set_xticks(range(len(properties)), [DISPLAY[p] for p in properties], rotation=25, ha="right")
    axs[1, 0].set_yticks(range(len(temps)), [f"{temp:.0f} K" for temp in temps])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            axs[1, 0].text(j, i, f"{matrix[i,j]:,}", ha="center", va="center", fontsize=8)
    axs[1, 0].set_title("(c) Common-temperature usable samples")
    fig.colorbar(im, ax=axs[1, 0], label="log10(N+1)")

    source_rows = [
        ("real ZT peak table", len(peak), "experiment; one row/sample"),
        ("sample relative density", metadata.relative_density_pct.notna().sum(), "experiment metadata"),
        ("sample grain size", metadata.grain_size_um.notna().sum(), "experiment metadata"),
        ("elastic-kL match", len(elastic), "cross-source formula match"),
        ("JARVIS 2D S/sigma", len(jarvis2d), "DFT; fixed n,T"),
    ]
    axs[1, 1].axis("off")
    axs[1, 1].text(0.02, 0.95, "(d) Additional data layers", transform=axs[1, 1].transAxes, va="top", weight="bold")
    y = 0.82
    for name, count, note in source_rows:
        axs[1, 1].text(0.04, y, f"{name}: {count:,}", transform=axs[1, 1].transAxes, weight="bold")
        axs[1, 1].text(0.04, y - 0.07, note, transform=axs[1, 1].transAxes, color="0.35")
        y -= 0.16
    axs[1, 1].text(
        0.04,
        0.02,
        "No direct wrinkle, pore geometry, phonon group velocity,\nor acoustic/optical-branch labels are present.",
        transform=axs[1, 1].transAxes,
        color="#8b1a1a",
        fontsize=8.5,
    )
    fig.suptitle("Data coverage: every relation uses its own maximum complete subset", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "01_data_coverage.png", bbox_inches="tight")
    plt.close(fig)
    return direct_counts, peak_pairs


def make_real_zt_figure(peak, manifest):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    panels = [
        ("carrier_concentration_cm3", "zt_peak", r"carrier concentration $n$ (cm$^{-3}$)", "ZT peak", True, False, "(a) Real ZT vs carrier concentration"),
        ("abs_seebeck_uV_K", "zt_peak", r"$|S|$ ($\mu$V/K)", "ZT peak", False, False, "(b) Real ZT vs Seebeck coefficient"),
        ("sigma", "zt_peak", r"electrical conductivity $\sigma$ (S/m)", "ZT peak", True, False, "(c) Real ZT vs electrical conductivity"),
        ("kappa_lattice", "zt_peak", r"lattice thermal conductivity $\kappa_L$ (W/mK)", "ZT peak", True, False, "(d) Real ZT vs lattice thermal conductivity"),
    ]
    for panel_id, (ax, spec) in enumerate(zip(axs.ravel(), panels), start=1):
        x, y, xlabel, ylabel, logx, logy, title = spec
        hb, n = density_panel(ax, peak, x, y, xlabel, ylabel, logx, logy, title)
        fig.colorbar(hb, ax=ax, label="log10 points/bin")
        manifest.append(
            {
                "figure": "02_experimental_ZT_global",
                "panel": chr(96 + panel_id),
                "relation": f"{y} vs {x}",
                "source": "Starrydata2 experimental curves",
                "n_observations": n,
                "n_unique_samples": n,
                "temperature_rule": "other property interpolated at each sample peak-ZT temperature",
                "true_ZT": True,
            }
        )
    fig.suptitle("All pairable experimental materials: one peak-ZT point per sample", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "02_experimental_ZT_global.png", bbox_inches="tight")
    plt.close(fig)


def temperature_scatter_panel(ax, frame, x, y, xlabel, ylabel, logx, logy, title, manifest, panel):
    colors = {300.0: "#2867b2", 600.0: "#d17a00", 900.0: "#b83280"}
    total = 0
    sample_ids = set()
    for temperature, group in frame.groupby("temperature_K"):
        subset = finite_pairs(group, x, y, positive_x=logx, positive_y=logy)
        total += len(subset)
        sample_ids.update(subset.sample_id.tolist())
        ax.scatter(
            subset[x],
            subset[y],
            s=9,
            alpha=0.28,
            edgecolors="none",
            color=colors.get(float(temperature), "grey"),
            label=f"{temperature:.0f} K (N={len(subset)})",
            rasterized=True,
        )
        add_binned_summary(ax, subset[x], subset[y], log_x=logx, bins=8, color=colors.get(float(temperature), "grey"))
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7)
    manifest.append(
        {
            "figure": "03_electronic_all_available",
            "panel": panel,
            "relation": f"{y} vs {x}",
            "source": "Starrydata2 experimental curves",
            "n_observations": total,
            "n_unique_samples": len(sample_ids),
            "temperature_rule": "same-sample interpolation at 300/600/900 K",
            "true_ZT": False,
        }
    )


def make_electronic_figure(common, manifest):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    specs = [
        ("carrier_concentration_cm3", "abs_seebeck_uV_K", r"$n$ (cm$^{-3}$)", r"$|S|$ ($\mu$V/K)", True, False, "(a) Pisarenko distribution: all n-S pairs"),
        ("carrier_concentration_cm3", "sigma", r"$n$ (cm$^{-3}$)", r"$\sigma$ (S/m)", True, True, "(b) Conductivity vs carrier concentration"),
        ("carrier_concentration_cm3", "mobility_cm2_Vs", r"$n$ (cm$^{-3}$)", r"mobility $\mu$ (cm$^2$/Vs)", True, True, "(c) Mobility vs carrier concentration"),
        ("carrier_concentration_cm3", "power_factor_used_mW_mK2", r"$n$ (cm$^{-3}$)", r"PF (mW/mK$^2$)", True, True, "(d) Power factor vs carrier concentration"),
    ]
    for panel, (ax, spec) in zip("abcd", zip(axs.ravel(), specs)):
        temperature_scatter_panel(ax, common, *spec, manifest, panel)
    fig.suptitle("Electronic transport: use every same-temperature pair; ZT not required", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "03_electronic_all_available.png", bbox_inches="tight")
    plt.close(fig)


def make_thermal_figure(kappa_grid, common, elastic, manifest):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    stats = kappa_grid.groupby("temperature_K").kappa_lattice.agg(
        n="count",
        p10=lambda x: x.quantile(0.10),
        p25=lambda x: x.quantile(0.25),
        median="median",
        p75=lambda x: x.quantile(0.75),
        p90=lambda x: x.quantile(0.90),
    )
    ax = axs[0, 0]
    ax.fill_between(stats.index, stats.p10, stats.p90, color="#9ecae1", alpha=0.35, label="10–90%")
    ax.fill_between(stats.index, stats.p25, stats.p75, color="#3182bd", alpha=0.35, label="25–75%")
    ax.plot(stats.index, stats["median"], color="#08519c", lw=2, label="median")
    ax.set_yscale("log")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"$\kappa_L$ (W/mK)")
    ax.set_title("(a) All experimental lattice-conductivity curves")
    ax.legend(frameon=False)
    manifest.append(
        {"figure":"04_thermal_all_available","panel":"a","relation":"kappa_lattice vs temperature percentiles","source":"Starrydata2 experimental curves","n_observations":len(kappa_grid),"n_unique_samples":kappa_grid.sample_id.nunique(),"temperature_rule":"sample median on 25 K grid","true_ZT":False}
    )

    ax = axs[0, 1]
    thermal_pair = common.dropna(subset=["kappa_lattice", "kappa_total"])
    for temperature, group in thermal_pair.groupby("temperature_K"):
        ax.scatter(group.kappa_total, group.kappa_lattice, s=8, alpha=0.25, label=f"{temperature:.0f} K (N={len(group)})", rasterized=True)
    line = np.geomspace(0.03, 100, 100)
    ax.plot(line, line, "k--", lw=1, label=r"$\kappa_L=\kappa_{total}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\kappa_{total}$ (W/mK)"); ax.set_ylabel(r"$\kappa_L$ (W/mK)")
    ax.set_title("(b) Lattice vs total thermal conductivity")
    ax.legend(frameon=False, fontsize=7)
    add_stats(ax, thermal_pair, "kappa_total", "kappa_lattice", "sample_id")
    manifest.append(
        {"figure":"04_thermal_all_available","panel":"b","relation":"kappa_lattice vs kappa_total","source":"Starrydata2 experimental curves","n_observations":len(thermal_pair),"n_unique_samples":thermal_pair.sample_id.nunique(),"temperature_rule":"same-sample interpolation at 300/600/900 K","true_ZT":False}
    )

    ax = axs[1, 0]
    wf = common.dropna(subset=["kappa_electronic", "sigma"]).query("kappa_electronic > 0 and sigma > 0")
    for temperature, group in wf.groupby("temperature_K"):
        ax.scatter(group.sigma, group.kappa_electronic, s=9, alpha=0.3, label=f"{temperature:.0f} K (N={len(group)})", rasterized=True)
        sigma_line = np.geomspace(max(group.sigma.min(), 1), group.sigma.max(), 80)
        ax.plot(sigma_line, 2.0e-8 * sigma_line * temperature, lw=1, alpha=0.7)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\sigma$ (S/m)"); ax.set_ylabel(r"$\kappa_e$ (W/mK)")
    ax.set_title("(c) Electronic heat conduction and Wiedemann-Franz trend")
    ax.legend(frameon=False, fontsize=7)
    add_stats(ax, wf, "sigma", "kappa_electronic", "sample_id")
    manifest.append(
        {"figure":"04_thermal_all_available","panel":"c","relation":"kappa_electronic vs sigma","source":"Starrydata2 experimental curves","n_observations":len(wf),"n_unique_samples":wf.sample_id.nunique(),"temperature_rule":"same-sample interpolation at 300/600/900 K","true_ZT":False}
    )

    ax = axs[1, 1]
    cross = elastic.dropna(subset=["bulk_sound_speed_proxy_m_s", "kL_300"])
    ax.scatter(cross.bulk_sound_speed_proxy_m_s, cross.kL_300, s=28, alpha=0.7, c=np.log10(cross.B_kv), cmap="plasma", edgecolors="none")
    ax.set_yscale("log")
    ax.set_xlabel(r"sound-speed proxy $\sqrt{B/\rho}$ (m/s)")
    ax.set_ylabel(r"experimental $\kappa_L$ at 300 K (W/mK)")
    ax.set_title("(d) Elastic/sound-speed proxy vs experimental kappaL")
    add_stats(ax, cross, "bulk_sound_speed_proxy_m_s", "kL_300")
    manifest.append(
        {"figure":"04_thermal_all_available","panel":"d","relation":"experimental kappa_lattice vs sqrt(B/rho)","source":"Starrydata2 kL + JARVIS 3D elastic crossmatch","n_observations":len(cross),"n_unique_samples":len(cross),"temperature_rule":"kL around 300 K; composition-level match","true_ZT":False}
    )
    fig.suptitle("Thermal transport: full experimental subsets plus explicit elastic crossmatch", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "04_thermal_all_available.png", bbox_inches="tight")
    plt.close(fig)


def boxplot_with_counts(ax, frame, category, value, order, title, ylabel):
    groups = [frame.loc[frame[category].eq(label), value].dropna().to_numpy() for label in order]
    ax.boxplot(groups, tick_labels=[f"{label}\nN={len(group)}" for label, group in zip(order, groups)], showfliers=False, patch_artist=True, boxprops=dict(facecolor="#9ecae1"), medianprops=dict(color="#8b1a1a", lw=1.5))
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=25)
    return sum(map(len, groups))


def make_structure_metadata_figure(peak, sample_meta, manifest):
    frame = peak.merge(sample_meta, on="sample_id", how="left")
    frame.to_csv(OUT / "experimental_ZT_with_structure_metadata.csv", index=False)
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.8))

    porosity = frame.dropna(subset=["porosity_fraction", "zt_peak"]).query("0 <= porosity_fraction <= 0.8").copy()
    bins = [-1e-9, 0.05, 0.10, 0.20, 0.40, 0.80]
    labels = ["0–5%", "5–10%", "10–20%", "20–40%", "40–80%"]
    porosity["porosity_bin"] = pd.cut(porosity.porosity_fraction, bins=bins, labels=labels)
    used = boxplot_with_counts(axs[0,0], porosity, "porosity_bin", "zt_peak", labels, "(a) Real peak ZT by porosity bin", "ZT peak")
    manifest.append({"figure":"05_structure_metadata_experimental","panel":"a","relation":"peak ZT by porosity bin","source":"Starrydata2 curves + sample metadata","n_observations":used,"n_unique_samples":used,"temperature_rule":"one peak-ZT value per sample","true_ZT":True})

    grain = frame.dropna(subset=["grain_size_um", "zt_peak"]).query("0.003 <= grain_size_um <= 500").copy()
    grain_bins = [0.003, 0.1, 1, 10, 100, 500]
    grain_labels = ["<0.1 um", "0.1–1 um", "1–10 um", "10–100 um", ">100 um"]
    grain["grain_bin"] = pd.cut(grain.grain_size_um, bins=grain_bins, labels=grain_labels, include_lowest=True)
    used = boxplot_with_counts(axs[0,1], grain, "grain_bin", "zt_peak", grain_labels, "(b) Real peak ZT by grain-size bin", "ZT peak")
    manifest.append({"figure":"05_structure_metadata_experimental","panel":"b","relation":"peak ZT by grain size","source":"Starrydata2 curves + sample metadata","n_observations":used,"n_unique_samples":used,"temperature_rule":"one peak-ZT value per sample","true_ZT":True})

    form = frame.dropna(subset=["form", "zt_peak"]).copy()
    counts = form.form.value_counts()
    order = counts[counts >= 30].head(8).index.tolist()
    used = boxplot_with_counts(axs[1,0], form, "form", "zt_peak", order, "(c) Real peak ZT by reported sample form", "ZT peak")
    manifest.append({"figure":"05_structure_metadata_experimental","panel":"c","relation":"peak ZT by sample form","source":"Starrydata2 curves + sample metadata","n_observations":used,"n_unique_samples":used,"temperature_rule":"one peak-ZT value per sample","true_ZT":True})

    density_k = frame.dropna(subset=["relative_density_pct", "kappa_lattice"]).query("20 <= relative_density_pct <= 100")
    hb = axs[1,1].hexbin(density_k.relative_density_pct, density_k.kappa_lattice, gridsize=30, mincnt=1, bins="log", yscale="log", cmap="viridis")
    add_binned_summary(axs[1,1], density_k.relative_density_pct, density_k.kappa_lattice, bins=8, color="white")
    axs[1,1].set_xlabel("reported relative density (%)")
    axs[1,1].set_ylabel(r"$\kappa_L$ (W/mK, log scale)")
    axs[1,1].set_title("(d) Relative density vs lattice thermal conductivity")
    add_stats(axs[1,1], density_k, "relative_density_pct", "kappa_lattice", "sample_id")
    fig.colorbar(hb, ax=axs[1,1], label="log10 points/bin")
    manifest.append({"figure":"05_structure_metadata_experimental","panel":"d","relation":"kappa_lattice vs relative density","source":"Starrydata2 curves + sample metadata","n_observations":len(density_k),"n_unique_samples":density_k.sample_id.nunique(),"temperature_rule":"kappaL interpolated at sample peak-ZT temperature","true_ZT":False})
    fig.suptitle("Experimental structure metadata: descriptive associations, not controlled causality", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "05_structure_metadata_experimental.png", bbox_inches="tight")
    plt.close(fig)


def jarvis_relation_panel(ax, frame, x, y, xlabel, ylabel, title, logy, manifest, panel):
    colors = {"n": "#2867b2", "p": "#c43d3d"}
    total = 0
    for carrier, group in frame.groupby("carrier"):
        subset = finite_pairs(group, x, y, positive_y=logy)
        total += len(subset)
        ax.scatter(subset[x], subset[y], s=8, alpha=0.25, color=colors[carrier], label=f"{carrier}-type N={len(subset)}", edgecolors="none", rasterized=True)
        add_binned_summary(ax, subset[x], subset[y], bins=10, color=colors[carrier])
    if logy: ax.set_yscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title); ax.legend(frameon=False, fontsize=7)
    manifest.append({"figure":"06_jarvis2d_shape_transport","panel":panel,"relation":f"{y} vs {x}","source":"JARVIS dft_2d","n_observations":total,"n_unique_samples":frame.dropna(subset=[x,y]).jid.nunique(),"temperature_rule":"fixed n=1e20 cm^-3 and T=600 K in source transport table","true_ZT":False})


def make_jarvis2d_figure(frame, manifest):
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 8.8))
    specs = [
        ("out_of_plane_ratio", "abs_seebeck_uV_K", "out-of-plane atomic span / sqrt(area)", r"$|S|$ ($\mu$V/K)", "(a) Atomic out-of-plane shape vs Seebeck", False),
        ("out_of_plane_ratio", "sigma_S_m", "out-of-plane atomic span / sqrt(area)", r"$\sigma$ (S/m)", "(b) Atomic out-of-plane shape vs conductivity", True),
        ("area_per_atom_A2", "abs_seebeck_uV_K", r"in-plane area per atom ($\AA^2$)", r"$|S|$ ($\mu$V/K)", "(c) In-plane packing proxy vs Seebeck", False),
        ("area_per_atom_A2", "sigma_S_m", r"in-plane area per atom ($\AA^2$)", r"$\sigma$ (S/m)", "(d) In-plane packing proxy vs conductivity", True),
    ]
    for panel, (ax, spec) in zip("abcd", zip(axs.ravel(), specs)):
        jarvis_relation_panel(ax, frame, *spec, manifest, panel)
    fig.suptitle("JARVIS 2D geometry vs transport: all available records; not real ZT", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "06_jarvis2d_shape_transport.png", bbox_inches="tight")
    plt.close(fig)


def write_summary(peak, common, metadata, elastic, jarvis2d, manifest, direct_counts, peak_pairs):
    lines = [
        "# 全体材料总体数据图：结果与边界",
        "",
        "## 核心修正",
        "",
        "原模型图只代表一个基准参数组合，不能解释为材料总体分布。本目录的主图现已改为：每个关系使用它所需的最大成对完整样本集；真实 ZT 面板每个实验样品只出现一次，其他物性在同一样品的峰值 ZT 温度配对。",
        "",
        "## 数据覆盖",
        "",
        "| 性质 | 直接实验样品数 |",
        "|---|---:|",
    ]
    for prop, count in sorted(direct_counts.items(), key=lambda item: -item[1]):
        lines.append(f"| {DISPLAY[prop]} | {count:,} |")
    lines += ["", "峰值 ZT 配对数：", "", "| 关系 | 同样品配对数 |", "|---|---:|"]
    for relation, count in peak_pairs.items():
        lines.append(f"| {relation} | {count:,} |")

    correlations = []
    for x, label in [
        ("carrier_concentration_cm3", "n"),
        ("abs_seebeck_uV_K", "|S|"),
        ("sigma", "sigma"),
        ("kappa_lattice", "kappaL"),
    ]:
        subset = peak.dropna(subset=[x, "zt_peak"])
        rho = spearmanr(subset[x], subset.zt_peak).statistic if len(subset) >= 3 else np.nan
        correlations.append((label, len(subset), rho))
    lines += [
        "",
        "## 真实 ZT 总体关系",
        "",
        "| 横轴 | N | Spearman rho(ZT,x) |",
        "|---|---:|---:|",
    ]
    for label, n, rho in correlations:
        lines.append(f"| {label} | {n:,} | {rho:+.3f} |")
    lines += [
        "",
        "散点密度和分箱中位数才是主要读图对象，不能用单条 SPB 曲线代替总体数据。由于样品来自不同化学家族、温度、制备状态和论文，相关性是描述性的，不是受控因果效应。",
        "电导配对优先使用直接电导曲线；样品只有电阻率曲线时使用 sigma=1/rho 补充，因此电导配对数可以高于直接电导覆盖数。",
        "",
        "## 结构数据能回答什么",
        "",
        f"- 相对密度可解析样品：{metadata.relative_density_pct.notna().sum():,}；可近似得到孔隙率，但许多记录是阈值或区间。",
        f"- 晶粒尺寸可解析样品：{metadata.grain_size_um.notna().sum():,}；用于分档，不把文本中值当高精度测量。",
        f"- 弹性—实验 κL 有效跨库配对：{len(elastic):,}；可检验 sqrt(B/rho) 声速代理，但存在化学式多晶型歧义。",
        f"- JARVIS 2D 几何—输运记录：{len(jarvis2d):,} 行；用于形状代理与 S/sigma，不是真实 ZT。",
        "- 当前数据没有显式褶皱幅度/波长、孔洞形貌、声子群速度或完整声学/光学支标签，因此不能把这些变量画成经验总体关系；模型图只作为待验证假设保留。",
        "",
        "## 图清单",
        "",
        "1. `01_data_coverage.png`：每种性质及每个配对关系到底用了多少样品；",
        "2. `02_experimental_ZT_global.png`：真实 ZT 与 n、S、sigma、kappaL 的全体可配对样品；",
        "3. `03_electronic_all_available.png`：不要求 ZT 的 Pisarenko、电导、迁移率、PF 总体数据；",
        "4. `04_thermal_all_available.png`：kappaL 温度带、总/晶格热导、电子热导及声速代理；",
        "5. `05_structure_metadata_experimental.png`：真实 ZT 与孔隙率、晶粒、样品形态；",
        "6. `06_jarvis2d_shape_transport.png`：二维原子几何代理与 S/sigma。",
        "",
        "每个面板的来源、样本数和温度规则见 `panel_manifest.csv`。",
    ]
    (OUT / "empirical_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    pd.DataFrame(manifest).to_csv(OUT / "panel_manifest.csv", index=False)


def main():
    print("loading and cleaning experimental curves...")
    curves, curve_metadata = load_curve_store()
    print("building peak-ZT and common-temperature pair tables...")
    peak = build_peak_table(curves, curve_metadata)
    common = build_common_temperature_table(curves, curve_metadata)
    kappa_grid = build_kappa_temperature_band(curves)
    sample_metadata = load_sample_metadata()
    elastic = load_elastic_kappa_crossmatch()
    jarvis2d = build_jarvis2d_geometry()

    manifest = []
    direct_counts, peak_pairs = make_coverage_figure(curves, peak, common, sample_metadata, elastic, jarvis2d)
    make_real_zt_figure(peak, manifest)
    make_electronic_figure(common, manifest)
    make_thermal_figure(kappa_grid, common, elastic, manifest)
    make_structure_metadata_figure(peak, sample_metadata, manifest)
    make_jarvis2d_figure(jarvis2d, manifest)
    write_summary(peak, common, sample_metadata, elastic, jarvis2d, manifest, direct_counts, peak_pairs)
    print(f"written: {OUT / 'empirical_summary.md'}")
    print(f"figures: {len(list(FIG.glob('*.png')))}; panels: {len(manifest)}")


if __name__ == "__main__":
    main()
