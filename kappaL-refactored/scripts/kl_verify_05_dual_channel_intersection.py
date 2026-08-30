"""Audit and intersect electronic power-factor and lattice-kappa channels.

Uses the existing 137 formula-matched starrydata2/JARVIS-3D table.  JARVIS PF
is available for 85 of those materials.  PF is maximized and experimental
kappa_L@300 K is minimized; no ZT is calculated because the temperatures and
data provenance differ.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core import Composition
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]


def canon(formula: str) -> str | None:
    try:
        return Composition(formula).reduced_formula
    except Exception:
        return None


def finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    return np.nan if value <= -99998 else value


def pareto_mask(pf: np.ndarray, kl: np.ndarray) -> np.ndarray:
    """Non-dominated points for max(PF), min(kL)."""
    keep = np.ones(len(pf), dtype=bool)
    for i, (pfi, kli) in enumerate(zip(pf, kl)):
        dominates = (pf >= pfi) & (kl <= kli) & ((pf > pfi) | (kl < kli))
        keep[i] = not dominates.any()
    return keep


def main() -> None:
    views = pd.read_parquet(ROOT / "features/kl_verify/kl_views.parquet")
    raw = json.load(open(ROOT / "data/raw/external/jarvis_kl/jdft_3d-8-18-2021.json"))
    by_jid = {row["jid"]: row for row in raw}

    by_formula = defaultdict(list)
    for row in raw:
        key = canon(row.get("formula", ""))
        if key:
            by_formula[key].append(row)

    ambiguity_rows = []
    for _, row in views.iterrows():
        polymorphs = by_formula[row.canon]
        energies = [finite(x.get("formation_energy_peratom")) for x in polymorphs]
        good = sorted(x for x in energies if np.isfinite(x))
        ambiguity_rows.append(
            {
                "canon": row.canon,
                "selected_jid": row.jid,
                "n_jarvis_polymorphs": len(polymorphs),
                "formation_energy_span_eV_atom": good[-1] - good[0] if len(good) > 1 else 0.0,
            }
        )
    ambiguity = pd.DataFrame(ambiguity_rows)
    ambiguity.to_csv(ROOT / "data/audit/kl_formula_polymorph_ambiguity.csv", index=False)

    all_rows = []
    summaries = []
    for carrier in ["n", "p"]:
        pf_field = f"{carrier}-powerfact"
        rows = []
        for _, row in views.iterrows():
            pf = finite(by_jid[row.jid].get(pf_field))
            if not np.isfinite(pf) or pf <= 0 or row.kL_300 <= 0:
                continue
            rows.append(
                {
                    "carrier": carrier,
                    "jid": row.jid,
                    "formula": row.formula,
                    "PF_jarvis": pf,
                    "kL_exp_300K": row.kL_300,
                    "Eg_opt": row.Eg_opt,
                    "B_kv": row.B_kv,
                    "G_gv": row.G_gv,
                }
            )
        frame = pd.DataFrame(rows)
        frame["logPF"] = np.log10(frame.PF_jarvis)
        frame["logkL"] = np.log10(frame.kL_exp_300K)
        frame["PF_percentile"] = frame.PF_jarvis.rank(pct=True)
        frame["low_kL_percentile"] = 1.0 - frame.kL_exp_300K.rank(pct=True) + 1.0 / len(frame)
        # Quantile thresholds give exactly the upper/lower fifth for N=85;
        # rank(pct)>=0.8 would include 18 rather than 17 observations.
        frame["top20_intersection"] = (
            frame.PF_jarvis >= frame.PF_jarvis.quantile(0.80)
        ) & (frame.kL_exp_300K <= frame.kL_exp_300K.quantile(0.20))
        frame["pareto"] = pareto_mask(
            frame.PF_jarvis.to_numpy(), frame.kL_exp_300K.to_numpy()
        )
        rho, pvalue = stats.spearmanr(frame.logPF, frame.logkL)
        summaries.append(
            {
                "carrier": carrier,
                "N": len(frame),
                "spearman_logPF_logkL": rho,
                "p": pvalue,
                "n_top20_intersection": int(frame.top20_intersection.sum()),
                "independent_expectation": len(frame) * 0.04,
                "n_pareto": int(frame.pareto.sum()),
            }
        )
        all_rows.append(frame)

    result = pd.concat(all_rows, ignore_index=True)
    result.to_csv(ROOT / "data/processed/pf_kL_dual_channel_intersection.csv", index=False)
    pd.DataFrame(summaries).to_csv(ROOT / "data/audit/pf_kL_dual_channel_summary.csv", index=False)
    result[result.pareto].to_csv(ROOT / "data/processed/pf_kL_pareto_front.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), constrained_layout=True)
    for ax, carrier in zip(axes, ["n", "p"]):
        frame = result[result.carrier == carrier]
        base = ax.scatter(
            frame.kL_exp_300K,
            frame.PF_jarvis,
            c=frame.Eg_opt,
            cmap="viridis",
            s=34,
            alpha=0.72,
            linewidths=0,
        )
        hits = frame[frame.top20_intersection]
        ax.scatter(
            hits.kL_exp_300K,
            hits.PF_jarvis,
            marker="*",
            s=150,
            facecolors="none",
            edgecolors="crimson",
            linewidths=1.3,
            label="top-PF / low-kL intersection",
        )
        front = frame[frame.pareto].sort_values("kL_exp_300K")
        ax.plot(front.kL_exp_300K, front.PF_jarvis, color="crimson", lw=1.0, label="Pareto front")
        for _, row in hits.iterrows():
            ax.annotate(row.formula, (row.kL_exp_300K, row.PF_jarvis), xytext=(4, 4), textcoords="offset points", fontsize=8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("experimental kappa_L at 300 K [W m-1 K-1]")
        ax.set_ylabel("JARVIS PF at 600 K, 1e20 cm-3 [database unit]")
        ax.set_title(f"{carrier}-type: electronic / lattice-channel intersection")
        ax.legend(frameon=False, fontsize=8)
        cb = fig.colorbar(base, ax=ax, pad=0.02)
        cb.set_label("OptB88vdW gap [eV]")
    fig.savefig(ROOT / "figures/pf_kL_dual_channel_intersection.png", dpi=220)
    fig.savefig(ROOT / "figures/pf_kL_dual_channel_intersection.pdf")
    plt.close(fig)

    print(pd.DataFrame(summaries).to_string(index=False))
    print("\nTop-20 intersections:\n", result[result.top20_intersection].to_string(index=False))
    print("\nPareto front:\n", result[result.pareto].to_string(index=False))
    print(
        f"\nFormula-polymorph ambiguity: {(ambiguity.n_jarvis_polymorphs > 1).sum()}/{len(ambiguity)} "
        f"matched formulas have >1 JARVIS structure; max={ambiguity.n_jarvis_polymorphs.max()}."
    )


if __name__ == "__main__":
    main()
