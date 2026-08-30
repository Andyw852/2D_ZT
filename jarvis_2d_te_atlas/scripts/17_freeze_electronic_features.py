"""L0-H: 冻结 Electronic View Features。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))
opt = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
opt_by_jid = {r["attributes"]["_jarvis_jid"]: r["attributes"] for r in opt}

def g(jid, f):
    a = opt_by_jid.get(jid, {}).get(f)
    return a if a not in (None, -99999, -99999.0) else np.nan

rows = []
for r in opt:
    jid = r["attributes"]["_jarvis_jid"]
    row = {"jid": jid}
    row["Eg_optb88vdw"] = g(jid, "_jarvis_optb88vdw_bandgap")
    row["Eg_mbj"] = g(jid, "_jarvis_mbj_bandgap")
    row["Eg_hse"] = g(jid, "_jarvis_hse_gap")
    # effective mass spectrum
    if jid in raw and "electron_mass_300K" in raw[jid]:
        em = np.sort(np.array([float(x) for x in raw[jid]["electron_mass_300K"]]))
        hm = np.sort(np.array([float(x) for x in raw[jid]["hole_mass_300K"]]))
        row["m_elec_median"] = np.median(em)
        row["m_elec_dom_geo"] = np.sqrt(em[0]*em[1]) if em[0] > 0 else np.nan
        row["m_elec_spectral_ratio"] = np.log10(em[2]/em[0]) if em[0] > 0 else np.nan
        row["m_hole_median"] = np.median(hm)
        row["m_hole_dom_geo"] = np.sqrt(hm[0]*hm[1]) if hm[0] > 0 else np.nan
        row["m_hole_spectral_ratio"] = np.log10(hm[2]/hm[0]) if hm[0] > 0 else np.nan
    else:
        for c in ["m_elec_median","m_elec_dom_geo","m_elec_spectral_ratio","m_hole_median","m_hole_dom_geo","m_hole_spectral_ratio"]:
            row[c] = np.nan
    rows.append(row)
df = pd.DataFrame(rows).sort_values("jid").reset_index(drop=True)
df.to_parquet(root / "features" / "electronic" / "electronic_features_v1.parquet", index=False)
print(f"electronic_features_v1: {df.shape}")
print(f"  Eg_optb88vdw coverage: {df['Eg_optb88vdw'].notna().sum()}/{len(df)}")
print(f"  m_elec_median coverage: {df['m_elec_median'].notna().sum()}/{len(df)}")
print(f"  Eg_mbj coverage: {df['Eg_mbj'].notna().sum()}, Eg_hse: {df['Eg_hse'].notna().sum()}")
print(f"  m_elec_median: median={df['m_elec_median'].median():.4f} p10={df['m_elec_median'].quantile(0.1):.4f} p90={df['m_elec_median'].quantile(0.9):.4f}")
print(f"  m_hole_median: median={df['m_hole_median'].median():.4f} p10={df['m_hole_median'].quantile(0.1):.4f} p90={df['m_hole_median'].quantile(0.9):.4f}")

meta = [
    ["Eg_optb88vdw", "OptB88vdW 带隙", "_jarvis_optb88vdw_bandgap", "raw", "eV", True, True, "100% 覆盖率，Electronic View 基础变量(0=金属)"],
    ["m_elec_median", "电子有效质量(稳健中位数, 面内类)", "electron_mass_300K 3 值", "median", "m_e", True, True, "mean 被面外大值污染~575x，median 才是物理值"],
    ["m_hole_median", "空穴有效质量(稳健中位数, 面内类)", "hole_mass_300K 3 值", "median", "m_e", True, True, "同上，mean 污染~300x"],
    ["m_elec_dom_geo", "电子有效质量两小主通道几何均值", "electron_mass_300K", "sqrt(s1*s2)", "m_e", True, False, "候选：面内类平均质量"],
    ["m_hole_dom_geo", "空穴有效质量两小主通道几何均值", "hole_mass_300K", "sqrt(s1*s2)", "m_e", True, False, "候选"],
    ["m_elec_spectral_ratio", "电子有效质量谱比值", "electron_mass_300K", "log10(max/min)", "log ratio", True, False, "候选：2D 质量各向异性"],
    ["m_hole_spectral_ratio", "空穴有效质量谱比值", "hole_mass_300K", "log10(max/min)", "log ratio", True, False, "候选"],
    ["Eg_mbj", "MBJ/TBmBJ 带隙", "_jarvis_mbj_bandgap", "raw", "eV", True, False, "22% 覆盖率，higher-level validation"],
    ["Eg_hse", "HSE06 带隙", "_jarvis_hse_gap", "raw", "eV", True, False, "4.9% 覆盖率，exploratory"],
]
mdf = pd.DataFrame(meta, columns=["feature_name","physical_meaning","source_property","transform","unit","permutation_invariant","used_in_main_model","reason"])
mdf.to_csv(root / "features" / "electronic" / "electronic_feature_metadata.csv", index=False)
print("\nwrote electronic features + metadata")
