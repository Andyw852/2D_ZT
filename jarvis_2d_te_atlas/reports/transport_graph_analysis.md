# Transport Graph 分析（Phase N）

## n/p Transport 图

- 严格使用 Phase L0 冻结的 6 个 V1 特征，RobustScaler + Euclidean，k=15。
- G_n_transport_v1：N=806，giant=1.0，单分量，0 isolated。
- G_p_transport_v1：N=803，giant=1.0，单分量，0 isolated。

## 关键修正

- 初版 D_sigma 在 s1=0（面外电导恰为 0）时产生 NaN，导致 16 个 n-type 节点孤立。
  改为 s1=0 时 D_sigma/A_sigma_dom 用 cap=6.0（高于有限值 max~5.5），表示"无限抑制"通道。修复后 0 isolated。

## Feature leave-one-out（n-type, k=15）

| 移除 feature | kNN overlap（越低影响越大） |
|---|---:|
| D_sigma | 0.596 |
| log_sigma_dom_geo | 0.677 |
| S_sign_fraction | 0.732 |
| A_sigma_dom | 0.726 |
| S_MAD | 0.744 |
| S_median | 0.787 |

- D_sigma（suppressed-channel contrast）与 log_sigma_dom_geo（主导通道尺度）对邻域影响最大。

## κ_e sensitivity

- σ-based 主图与 κ_e-based sensitivity 图（T2）距离相关 >0.96，kNN overlap ~0.72，
  确认 κ_e 只是 sensitivity，不进主模型。
