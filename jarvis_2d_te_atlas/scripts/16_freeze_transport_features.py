"""L0-G: 冻结最终 n/p Transport Features V1。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

def build(carrier):
    n = "n" if carrier == "n" else "p"
    sf, cf, kf, pf = n + "seeb", n + "cond", n + "kappa", n + "pf"
    rows = []
    for jid in raw:
        if not all(k in raw[jid] for k in (sf, cf, kf, pf)):
            continue
        S = np.array([float(x) for x in raw[jid][sf]])
        C = np.sort(np.array([float(x) for x in raw[jid][cf]]))
        K = np.sort(np.array([float(x) for x in raw[jid][kf]]))
        PF = np.array([float(x) for x in raw[jid][pf]])
        med = np.median(S); mad = np.median(np.abs(S - med))
        pos = int((S > 0).sum()) if carrier == "p" else int((S < 0).sum())
        s1, s2, s3 = C; k1, k2, k3 = K
        r = {"jid": jid}
        # Seebeck
        r["S_mean"] = S.mean(); r["S_median"] = med; r["S_std"] = S.std()
        r["S_MAD"] = mad; r["S_min"] = S.min(); r["S_max"] = S.max(); r["S_range"] = S.max() - S.min()
        r["S_abs_mean"] = np.abs(S).mean(); r["S_abs_max"] = np.abs(S).max()
        r["S_relative_spread"] = S.std() / (abs(S.mean()) + 1e-12)
        r["S_sign_fraction"] = pos / 3
        # sigma spectrum
        r["sigma_mean"] = C.mean(); r["sigma_median"] = np.median(C)
        r["sigma_dom_geo"] = np.sqrt(s2 * s3) if s2 > 0 else np.nan
        r["sigma_geo_all"] = (s1*s2*s3)**(1/3) if s1 > 0 else np.nan
        r["log_sigma_mean"] = np.log10(C.mean() + 1e-6)
        r["log_sigma_dom_geo"] = np.log10(np.sqrt(s2*s3) + 1e-6) if s2 > 0 else 0.0
        # 当 s1=0（面外电导恰为 0，即"无限抑制"）时，D/A 用大值 cap=6.0（高于有限值 max ~5.5）
        r["D_sigma"] = np.log10(s2/s1) if s1 > 1e-9 else 6.0
        r["A_sigma_dom"] = np.log10(s3/s2) if s2 > 1e-9 else 6.0
        r["A_sigma_total"] = np.log10(s3/s1) if s1 > 1e-9 else 6.0
        # kappa spectrum (sensitivity)
        r["kappa_dom_geo"] = np.sqrt(k2*k3) if k2 > 0 else np.nan
        r["log_kappa_dom_geo"] = np.log10(np.sqrt(k2*k3) + 1e-6) if k2 > 0 else 0.0
        r["D_kappa"] = np.log10(k2/k1) if k1 > 1e-9 else 6.0
        r["A_kappa_dom"] = np.log10(k3/k2) if k2 > 1e-9 else 6.0
        # PF (external label)
        r["PF_mean"] = PF.mean()
        rows.append(r)
    return pd.DataFrame(rows).sort_values("jid").reset_index(drop=True)

# 冻结的 V1 主模型特征
V1 = ["S_median", "S_MAD", "S_sign_fraction", "log_sigma_dom_geo", "D_sigma", "A_sigma_dom"]

metadata_rows = []
for carrier in ["n", "p"]:
    df = build(carrier)
    # v1 文件（jid + V1 主特征 + PF 外部标签）
    v1 = df[["jid"] + V1 + ["PF_mean"]]
    v1.to_parquet(root / "features" / "transport" / f"{carrier}_transport_features_v1.parquet", index=False)
    # candidates 文件（全部候选）
    df.to_parquet(root / "features" / "transport" / f"{carrier}_transport_features_candidates.parquet", index=False)
    print(f"{carrier}-type: v1={v1.shape} candidates={df.shape}")

# metadata
meta = [
    ["S_median", "Seebeck 主值中位数（稳健载流子符号/量级）", "n/p-Seebeck 3 本征值", "median", "uV/K", True, True, "对单个符号异常本征值稳健(100% 匹配多数符号)"],
    ["S_MAD", "Seebeck 主值谱稳健离散度", "n/p-Seebeck 3 本征值", "MAD", "uV/K", True, True, "稳健 spread，不被单个异常值主导"],
    ["S_sign_fraction", "Seebeck 主值符号一致比例(0-1)", "n/p-Seebeck 3 本征值", "count/3", "dimensionless", True, True, "捕获 bipolar/small-gap 符号混合(与 Eg 强相关)"],
    ["log_sigma_dom_geo", "两个较强电导主通道几何尺度(对数)", "sigma 3 本征值", "log10(sqrt(s2*s3))", "log(S/m)", True, True, "主导输运尺度，避免 weakest channel 污染"],
    ["D_sigma", "suppressed-channel contrast", "sigma 3 本征值", "log10(s2/s1)", "log ratio", True, True, "quasi-2D 维度信息(主通道/被抑制通道)"],
    ["A_sigma_dom", "dominant-channel anisotropy", "sigma 3 本征值", "log10(s3/s2)", "log ratio", True, True, "两主导通道差异(多数材料≈0)"],
    ["log_kappa_dom_geo", "kappa_e 主导通道尺度(对数)", "kappa_e 3 本征值", "log10(sqrt(k2*k3))", "log ratio", True, False, "与 sigma 高度冗余，作 sensitivity"],
    ["D_kappa", "kappa_e suppressed-channel contrast", "kappa_e 3 本征值", "log10(k2/k1)", "log ratio", True, False, "sensitivity"],
    ["A_kappa_dom", "kappa_e dominant-channel anisotropy", "kappa_e 3 本征值", "log10(k3/k2)", "log ratio", True, False, "sensitivity"],
    ["PF_mean", "功率因子(数据库定义)", "PF 3 本征值", "mean", "DB unit", True, False, "external performance label, 不进 embedding"],
    ["S_mean", "Seebeck 算术均值", "Seebeck 3 本征值", "mean", "uV/K", True, False, "被符号异常本征值污染，不用于主模型"],
]
mdf = pd.DataFrame(meta, columns=["feature_name","physical_meaning","source_property","transform","unit","permutation_invariant","used_in_main_model","reason"])
for carrier in ["n", "p"]:
    mdf.to_csv(root / "features" / "transport" / f"{carrier}_transport_feature_metadata.csv", index=False)
print("\nwrote v1 parquet + candidates + metadata")
print("Frozen V1 features:", V1)
