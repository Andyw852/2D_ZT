# Transport Representation 选择审计（L0-D / L0-F）

## R1 / R2 / R3 比较（L0-D）

- R1（mean-only）：[S_mean, log_sigma_mean]
- R2（full spectrum）：[S_median, S_MAD, S_range, log_sigma_mean, D_sigma, A_sigma_dom]
- R3（dominant-channel）：[S_median, S_MAD, log_sigma_dom_geo, A_sigma_dom, D_sigma]

| pair | n dist_spearman | n knn10 | p dist_spearman | p knn10 |
|---|---:|---:|---:|---:|
| R1 vs R2 | 0.528 | 0.152 | 0.586 | 0.139 |
| R1 vs R3 | 0.549 | 0.160 | 0.594 | 0.158 |
| R2 vs R3 | 0.920 | 0.547 | 0.945 | 0.567 |

结论：mean-only（R1）与 spectrum 表示（R2/R3）的邻域**差异巨大**（kNN 重叠仅 ~14–16%），
证实 OPTIMADE 标量均值丢失关键信息。R2 与 R3 较一致（dist 相关 0.92–0.94，kNN 重叠 ~0.55），
但仍非 >0.8，说明 weakest channel 的处理方式（R2 用 mean 含它，R3 用 dom_geo 去掉它）对分类有中等影响。

## T1 / T2 / T3 比较（L0-F）

- T1（主模型, sigma）：[S_median, S_MAD, S_range, log_sigma_dom_geo, D_sigma, A_sigma_dom]
- T2（kappa_e 替代）：同结构，sigma→kappa_e
- T3（sigma+kappa_e 全部）：T1 + kappa_e 谱量

| pair | n dist_spearman | n knn10 | p dist_spearman | p knn10 |
|---|---:|---:|---:|---:|
| T1 vs T2 | 0.972 | 0.721 | 0.971 | 0.715 |
| T1 vs T3 | 0.965 | 0.767 | 0.967 | 0.767 |
| T2 vs T3 | 0.975 | 0.797 | 0.979 | 0.787 |

| PCA effective rank | T1 | T2 | T3 |
|---|---:|---:|---:|
| n | 4.15 | 4.35 | 4.04 |
| p | 3.98 | 4.20 | 3.98 |

结论：
1. T1 / T2 / T3 的距离结构高度相关（Spearman 0.96–0.98），kNN 重叠 0.72–0.80，
   即 **使用 sigma 或 kappa_e 得到几乎相同的材料空间**。
2. 加入 kappa_e（T3）对 T1 的改变有限（kNN 重叠 0.77，未达 >0.9），属中等冗余。
3. 遵循「不重复计权同一输运自由度」原则，**主模型选 T1（sigma 版，更简洁），kappa_e 作 sensitivity（T2/T3）**。
4. Transport View effective dimension ≈ 4（T1 的 PCA effective rank 3.98–4.15），
   显著高于 mean-only 的 ~1.4，说明 spectrum 确实提供 ~3 个额外自由度。

## 最终选择原则（第 36 节）

选中 T1/R3 风格：permutation-invariant、对单异常主值稳健（median/MAD）、保留显著 anisotropy（D_sigma 为主）、
不被 weakest channel 完全支配（用 log_sigma_dom_geo 而非 mean）、与 PF landscape 关联、n/p 对称、特征简洁。
