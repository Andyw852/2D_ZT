# JARVIS 二维热电输运图谱 第二阶段（Phase G–K）执行报告

## 一、Phase G：原始 3 本征值恢复 —— 成功

**关键发现**：figshare 在本环境仍 403，但 JARVIS 官方 NIST 静态 XML 服务器可达：

    https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/{JID}

（OPTIMADE 每条记录的 _jarvis_reference 字段正是指向这里）。每个 JID 返回完整 XML，其中：

- <boltztrap_info> 含 p/n 的 seeb / cond / pf / kappa 各 **3 个本征值**；
- <effective_mass> 含 electron/hole_mass_300K 各 3 个本征值。

下载全部 1103 个 JID（0 失败）：**807 条含输运本征值，678 条含有效质量**，与 OPTIMADE coverage 一致。
复数检查：0 个复数 / 0 个非数值（全部实数）。

## 二、12 个必须回答的问题

1. **是否成功恢复原始 3 eigenvalues？** 是（807 条输运 + 678 条有效质量，全部 3 值实数）。
2. **是否能够获取原始 3×3 tensor？** 否。XML 只存已对角化的 3 本征值，不含 3×3 张量与本征向量，
   因此只能研究 principal-value spectrum，不能研究晶向输运方向（第 38/39 节）。
3. **OPTIMADE scalar 是否严格等于 3 eigenvalues 的 arithmetic mean？** **是，严格相等**
   （400 个随机样本 MAE=0、最大误差=0、400/400 精确匹配，见 reports/optimade_mean_validation.csv）。
4. **n-S、n-sigma、n-kappa_e 是否同一批 JID？** 是，**完全相同**：|S_n ∩ sigma_n ∩ kappa_e_n| = 806。
5. **p 型对应集合是否相同？** 是，p 型完全相同 803 条（n=806 vs p=803，相差 4 条 n-only、1 条 p-only）。
6. **sigma 与 kappa_e 相关系数？** n：Pearson 0.960 / Spearman 0.928；p：Pearson 0.965 / Spearman 0.947 —— 高度相关。
7. **PF 与基础输运量的冗余程度？** PF = S²σ/1e6（源码已确认）。与 S、与 σ 的单独相关性都弱（Pearson 0.05/0.19），
   不可由单一基础量替代；但它是派生量，应作 external label。
8. **Transport View 的 effective dimension？** mean-only ≈ 1.3–1.4（近 1 维，因 σ 与 κ_e 共线）；
   加入 tensor spectrum 后 ≈ 3.1–3.2（真实 3 维）。
9. **tensor spectral anisotropy 是否明显改变材料邻域？** **是，显著改变**：kNN(10) overlap 仅 18–20%，
   距离排序 Spearman 仅 0.62–0.66。anisotropy 是重要自由度，OPTIMADE 标量均值确实丢失了关键信息。
10. **最终推荐使用哪些 transport variables？** S_mean、log_sigma_mean、S_std、S_range、
    sigma_anisotropy_log、（可选 kappa_e_anisotropy_log）；PF_mean 仅作外部标签。
11. **最终推荐建立几个 physical views？** 4 个：Structure、Electronic、n-Transport、p-Transport。
12. **PF 是否继续作为 external performance label？** 是。

## 三、最终决策表（Phase K）

| Candidate View | Features | N | Keep? | Reason |
|---|---|---:|---|---|
| Structure | SOAP + composition + local geometry | 1103 | YES | backbone（最高覆盖率） |
| Electronic | OptB88vdW Eg + effective mass | 1103 (Eg) / 678 (mass) | YES | 电子结构层，mass 子集 enrich |
| n-Transport | S_mean + log_sigma_mean + spectrum | 806 | YES | 主输运层，保留 anisotropy |
| n-Transport extended | + log_kappa_e_mean | 806 | 可选(sensitivity) | sigma 与 kappa_e 高度冗余 |
| p-Transport | S_mean + log_sigma_mean + spectrum | 803 | YES | 主输运层 |
| PF-n | PF_n | 806 | NO as input | external validation label |
| PF-p | PF_p | 803 | NO as input | external validation label |
| kappa_L | unavailable | 0 | NO | not available |
| ZT | unavailable | 0 | NO | not available |

**冗余结论**：sigma 与 kappa_e 高度共线（r≈0.93–0.96），主模型默认只保留其一（sigma），
kappa_e 作为 Model T2 sensitivity；anisotropy 必须保留（否则丢失约 1 个有效维度、邻域错乱 80%）。

## 四、本轮生成的数据文件

- data/processed/transport_eigenvalues_raw.json（807 输运 + 678 有效质量，原始 3 本征值）
- data/processed/transport_eigenvalues.parquet（1103×31）
- features/transport/n_transport_tensor_features.parquet（806×26）
- features/transport/p_transport_tensor_features.parquet（803×26）
- data/audit/view_jid_overlap.csv
- data/audit/transport_numeric_audit.csv
- data/audit/transport_correlation_pearson.csv
- data/audit/transport_correlation_spearman.csv
- data/audit/transport_pca_diagnostics.csv
- data/audit/seebeck_sign_violations.csv
- reports/raw_transport_source_audit.md
- reports/jarvis_transport_source_probe.csv
- reports/optimade_mean_validation.csv
- reports/transport_redundancy_analysis.md
- reports/round2_summary.md（本文件）

## 五、STOP

按第 49 节要求，Phase K 完成后停止，未进行 SOAP / kNN graph / UMAP / Diffusion Map / manifold alignment /
supra graph / joint embedding / unified atlas。下一阶段（Phase L 起）需先据此确认 physical views 与输运表示。
