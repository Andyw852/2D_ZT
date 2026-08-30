# Phase L0 执行报告：二维输运主值谱物理一致性审计与最终特征冻结

## Executive Summary

本阶段从 JARVIS 的 3 个 tensor principal values（无本征向量）中，提取了最少但物理意义最清晰的输运自由度，
并冻结了 V1 特征集。三大决定性结论：

1. **PF 的 eigenvalue pairing 高度敏感**（中位歧义 0.47–0.64），PF 只能作 database-defined external label，不进 embedding。
2. **Seebeck 用 median 而非 mean**（median 100% 匹配多数符号；mean 被 15% 的符号异常本征值污染）；
   符号异常强烈关联 small/zero band gap（p~1e-47），是 bipolar 输运而非数值噪声。
3. **conductivity 谱的 quasi-2D 特征明确**：suppressed-channel contrast D_sigma 很大（中位 ~3.1–3.3），
   dominant-channel anisotropy A_sigma_dom≈0；kappa_e 谱与 sigma 谱几乎一致（r~0.96–0.98），故 kappa_e 仅作 sensitivity。

## 17 个问题的回答

1. **PF 对 S 与 σ 的 eigenvalue pairing 是否敏感？** 是，高度敏感（中位歧义 n=0.637 / p=0.473，≥0.5 占 60%/48%）。
2. **PF 可否继续作可靠的 database 内部 external label？** 可以，但必须定义为 JARVIS pairing convention 下的 database-defined PF，不解释方向性。
3. **conductivity spectrum 是否普遍 weakest << two dominant？** 是，D_sigma 中位 3.33(n)/3.06(p)，即 weakest 比次弱小 ~1000–2000 倍。
4. **full-spectrum anisotropy 主要反映哪个？** 主要反映 suppressed-channel contrast (D_sigma)；A_sigma_dom 中位≈0。
5. **两个 dominant 通道差异多大？** 很小，A_sigma_dom 中位 0（p90 仅 0.74，~5.5 倍）。
6. **κ_e 与 σ 的 spectrum anisotropy 是否高度一致？** 是，D 与 A_dom 的 r≈0.95–0.98。
7. **15% Seebeck 符号违例是否与 small gap 显著相关？** 是，bad-sign 材料 gap 中位 0.000 eV vs good-sign 1.34–1.36 eV，Mann-Whitney p~1e-47。
8. **S_median 是否比 S_mean 更符合多数符号？** 是，median 100% vs mean 92.8–93.2%。
9. **R1/R2/R3 邻域差别？** R1 vs R2/R3 knn10≈0.14–0.16（巨大差异）；R2 vs R3 knn10≈0.55。
10. **anisotropy spectrum 信息是否应保留？** 是，必须保留（否则丢失主要输运自由度）。
11. **T1/T2/T3 哪个作主模型？** T1（sigma 版，简洁且与 T2/T3 几乎等价）。
12. **κ_e 是否只作 sensitivity？** 是（T1 vs T2 knn10 0.72、dist 相关 0.97）。
13. **PF 是否继续不进主 Transport View？** 是，external label。
14. **effective-mass spectrum 如何表示？** median（稳健面内类）为主，dom_geo/spectral_ratio 为候选；不用 mean（被污染 300–575 倍）。
15. **最终 Electronic View 变量？** Eg_optb88vdw + m_elec_median + m_hole_median（rich 层 678；Eg 层 1103）。
16. **最终 n-Transport View 变量？** S_median, S_MAD, S_sign_fraction, log_sigma_dom_geo, D_sigma, A_sigma_dom。
17. **最终 p-Transport View 变量？** 与 n 对称：S_median, S_MAD, S_sign_fraction, log_sigma_dom_geo, D_sigma, A_sigma_dom。

## 最终 Feature Decision Table

| View | Feature | Keep in main model? | Role | Reason |
|---|---|---|---|---|
| n-Transport | S_median | KEEP | primary | 稳健符号/量级，100% 匹配多数符号 |
| n-Transport | S_MAD | KEEP | dispersion | 稳健 spread |
| n-Transport | S_sign_fraction | KEEP | sign cleanliness | bipolar/small-gap 标志，与 Eg 强相关但非 Eg |
| n-Transport | log_sigma_dom_geo | KEEP | transport scale | 主导通道尺度，不被 weakest 污染 |
| n-Transport | D_sigma | KEEP | dimensionality contrast | quasi-2D 特征 |
| n-Transport | A_sigma_dom | KEEP | dominant anisotropy | 主导通道差异（多数≈0） |
| n-Transport | log_kappa_dom_geo / D_kappa / A_kappa_dom | SENSITIVITY | validation | 与 sigma 冗余(r~0.96-0.98) |
| n-Transport | PF_mean | EXTERNAL LABEL | performance | 派生量 + pairing 敏感，不进 embedding |
| n-Transport | S_mean | DROP | — | 被符号异常本征值污染 |
| p-Transport | (与 n 对称同 6 个) | KEEP | 同上 | 同上 |
| Electronic | Eg_optb88vdw | KEEP | primary | 100% 覆盖，0=金属 |
| Electronic | m_elec_median | KEEP | primary | 稳健面内类有效质量 |
| Electronic | m_hole_median | KEEP | primary | 同上 |
| Electronic | m_*_dom_geo / m_*_spectral_ratio | CANDIDATE | spectrum | 2D 质量各向异性 |
| Electronic | Eg_mbj | VALIDATION | higher-level | 22% 覆盖 |
| Electronic | Eg_hse | EXPLORATORY | — | 4.9% 覆盖 |

## 冻结的 V1（FINAL_TRANSPORT_REPRESENTATION_V1）

n-Transport V1（806 条）= [S_median, S_MAD, S_sign_fraction, log_sigma_dom_geo, D_sigma, A_sigma_dom]
p-Transport V1（803 条）= [S_median, S_MAD, S_sign_fraction, log_sigma_dom_geo, D_sigma, A_sigma_dom]
Electronic V1（1103 条）= [Eg_optb88vdw, m_elec_median, m_hole_median]（m_* 覆盖 678）
PF = external performance label；κ_e = sensitivity；MBJ/HSE = validation/exploratory。

## 语言规范遵守（第 62 节）

全程使用 principal value / dominant channel / weakest channel / suppressed-channel contrast / dominant-channel anisotropy /
quasi-2D-like spectral signature，未使用 x/y/z 或 in-plane/out-of-plane 方向断言（无本征向量）。

## STOP

Phase L0 完成。未执行 SOAP / kNN graph / UMAP / t-SNE / Diffusion Map / Spectral Embedding / SNF / multilayer graph /
supra adjacency / manifold alignment / joint Laplacian / superlattice generation（第 66 节）。
