# 多视图 Layer Strength Normalization（P-B）

各 layer 缩放到 mean node strength ≈ 1（W_scaled = W / mean_strength），避免高密度 layer（如金属 clique 的 Eg 层）在 supra 图中过度主导。

| view | N | mean_strength_before | scale_factor |
|---|---:|---:|---:|
| structure_v1 | 1103 | 8.83 | 0.1132 |
| Eg_v1 | 1103 | 17.69 | 0.0565 |
| m_electron_v1 | 678 | 16.86 | 0.0593 |
| m_hole_v1 | 678 | 16.74 | 0.0597 |
| n_transport_v1 | 806 | 9.18 | 0.1089 |
| p_transport_v1 | 803 | 9.17 | 0.1090 |

- Eg 与 mass 层 mean strength 更高（金属 clique 中 Eg=0 → affinity=1 产生高密度连接）。
- alpha 全部 = 1（第一版），未用 PF 调权。
