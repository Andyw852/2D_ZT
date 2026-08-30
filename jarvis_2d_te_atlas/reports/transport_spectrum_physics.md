# 输运主值谱物理一致性审计（L0-B / L0-E）

## Conductivity 谱的二维性

将 3 个本征值按数值排序 sigma_(1) <= sigma_(2) <= sigma_(3)：

| 量 | n-type | p-type |
|---|---:|---:|
| D_sigma = log10(s2/s1) 中位数 | 3.333 (p10 2.197, p90 4.804) | 3.061 (p10 2.013, p90 4.475) |
| A_sigma_dom = log10(s3/s2) 中位数 | 0.000 (p90 0.737) | 0.000 (p90 0.817) |
| A_sigma_total = log10(s3/s1) 中位数 | 3.521 | 3.309 |

结论：
1. **D_sigma 很大（中位 ~3.1–3.3，即 weakest channel 比次弱小 ~1000–2000 倍）**，
   普遍存在 suppressed channel，即 quasi-2D-like 谱特征明显。
2. **A_sigma_dom 很小（中位 0）**，两个 dominant channels 近乎简并，
   面内类输运接近各向同性（p90 仅 ~0.74，即 ~5.5 倍）。
3. 因此 full-spectrum anisotropy 主要由 **suppressed-channel contrast (D_sigma)** 主导，
   而非 dominant-channel anisotropy (A_sigma_dom)。

## kappa_e 谱与 sigma 谱的一致性（第 16 节）

| 量 | n Pearson | n Spearman | p Pearson | p Spearman |
|---|---:|---:|---:|---:|
| D (suppressed contrast) | 0.965 | 0.965 | 0.957 | 0.955 |
| A_dom (dominant anisotropy) | 0.975 | 0.979 | 0.971 | 0.981 |

结论：kappa_e 与 sigma 的谱形状（D 与 A_dom）高度一致（r≈0.96–0.98），
进一步证明 kappa_e 是 sigma 的冗余 sensitivity variable，不是独立主自由度。

## Effective-mass 谱（L0-E）

| 量 | electron (n=678) | hole (n=678) |
|---|---:|---:|
| mean 中位数 | 71.59（p90 10113） | 110.53（p90 15749） |
| median 中位数 | **0.560**（p90 5.46） | **1.135**（p90 10.77） |
| spectral_ratio 中位数 | 3.207 | 3.069 |
| 负本征值 | 0 | 0 |
| median(mean/median) | 574.6 | 301.2 |

结论：effective mass 的 arithmetic mean 被面外类大值（~1e3–1e4）严重污染（mean/median ≈ 300–575 倍），
OPTIMADE 的 avg_elec/avg_hole_mass 因此不是物理有效质量。**median（electron 0.56 m_e、hole 1.14 m_e）才是
稳健的面内类有效质量**。全部本征值为正（无需 abs）。
