# 输运冗余性与 anisotropy 审计报告（Phase J）

## 1. sigma 与 kappa_e 的相关性（第 32 节核心问题）

| carrier | Pearson(log_sigma, log_kappa_e) | Spearman(log_sigma, log_kappa_e) |
|---|---|---|
| n | 0.9598 | 0.9278 |
| p | 0.9649 | 0.9474 |

结论：两者高度相关（Pearson > 0.95；Spearman 0.93–0.95）。物理原因是固定 T、固定 doping、常数 τ 下
电子热导与电导率近似满足 Wiedemann–Franz 关系 κ_e ≈ L·σ·T。**sigma 与 kappa_e 不应作为两个完全独立、
等权的 property layers**（第 32/42 节）。

## 2. PF 的派生冗余（第 33 节）

- 源码公式（逐 principal value）：PF_i = S_i² · σ_i / 1e6。用恢复的 3 本征值验证：
  4827 个本征值三元组，最大相对误差约 0.03%（仅来自 S、σ 入库前 round(k,2) 的舍入），公式**确认**。
- 相关性：
  - corr(PF, S)：n Pearson 0.053 / Spearman -0.550；p Pearson -0.025 / Spearman 0.445 —— 弱到中等。
  - corr(PF, log_sigma)：n/p Pearson 约 0.19 —— 弱。
- 含义：PF 由 S² 与 σ 的乘积决定，与二者单独都不强相关，**不能由 S 或 σ 单独替代**；
  但它是派生量，**默认作为 external performance label，不进入 embedding**（第 34 节）。

## 3. Transport View 的 effective dimension（第 36 节 PCA 诊断）

| representation | 特征 | n | PC1 | PC2 | PC3 | 累计(2维) | effective rank |
|---|---|---|---:|---:|---:|---:|---:|
| n mean-only | [S_mean, log_sigma, log_kappa_e] | 806 | 0.843 | 0.153 | 0.004 | 0.996 | 1.36 |
| p mean-only | 同上 | 803 | 0.863 | 0.133 | 0.004 | 0.996 | 1.31 |
| n spectrum | 7 个（含 std/range/anisotropy） | 790 | 0.422 | 0.291 | 0.216 | 0.713 | 3.19 |
| p spectrum | 同上 | 796 | 0.423 | 0.322 | 0.198 | 0.744 | 3.09 |

结论：mean-only 输运空间接近 **1 维**（effective rank 1.3–1.4，因 sigma/kappa_e 高度共线，
本质是 [S, 一个输运量] 二维，再加一点残差）。加入 tensor spectrum 后 effective rank 升到约 **3**，
说明 anisotropy/principal-value dispersion 是真实存在的第三自由度。

## 4. tensor spectral anisotropy 是否改变材料邻域（第 37 节）

| carrier | n | Spearman(dist_mean, dist_spectrum) | kNN(10) overlap |
|---|---:|---:|---:|
| n | 790 | 0.6245 | 0.201 |
| p | 796 | 0.6642 | 0.177 |

结论：**mean-only 与 tensor-spectrum 表示给出的材料邻域显著不同**——kNN(10) 重叠只有约 18–20%，
距离排序相关系数约 0.62–0.66。即 **anisotropy 是重要输运自由度**，OPTIMADE 的标量均值确实丢失了
关键信息（第 3 节担心的数据损失成立）。

## 5. 数值质量与符号检查（Phase H2）

- sigma / kappa_e / PF：全部非负（sigma 有 16 个本征值恰为 0，kappa_e 1 个，PF 有 154 个为 0），无负值、无 NaN/Inf。
- Seebeck 符号：n 型应负、p 型应正。按 3 本征值均值判据，n 型 121/806（15.0%）、p 型 123/803（15.3%）出现
  符号违例；按“任一本征值违例”计 n=242、p=202。**未删除**，JID 已输出到 data/audit/seebeck_sign_violations.csv 供后续人工检查。
  初步判断这些违例主要来自面外（out-of-plane）本征值量级小、符号不稳定（二维材料的真空方向输运无物理意义）。

## 6. 最终推荐 transport variables（Phase K 依据）

1. 基础量：S_mean、log_sigma_mean（sigma 与 kappa_e 取其一，默认 sigma；kappa_e 作为 sensitivity Model T2）。
2. spectrum 量：S_std、S_range、sigma_anisotropy_log、kappa_e_anisotropy_log（anisotropy 显著改变邻域，必须保留）。
3. PF_mean 仅作 external label，不进 embedding。
4. 距离度量使用 log10 后的 sigma/kappa_e（跨数量级），epsilon 已按 min_positive/10 确定并记录
   （eps_sigma=0.0933，eps_kappa=4.643e9）。
