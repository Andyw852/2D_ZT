# External Label Diagnostics（Phase O）

> **注意**：PF_mean = S²σ 与 transport 图特征（S、σ 谱）同源，故 PF smoothness 是**同源一致性**，非独立外部验证。

## PF smoothness（log10(PF+eps)）+ null test（1000 次打乱）

| graph | smoothness | null mean±std | z | p |
|---|---:|---:|---:|---:|
| n_transport_v1 | 0.530 | 1.211±0.028 | -24.0 | 0.000 |
| p_transport_v1 | 0.541 | 1.176±0.025 | -25.8 | 0.000 |
| Structure（n 共同 JID） | 1.063 | — | -3.6 | 0.000 |
| Electronic_n | 0.994 | — | -1.0 | 0.089 |

- **PF 在 Transport graph 上高度平滑**（z≈-24～-26，远低于随机），但这是**同源一致性**（PF 与 transport 特征同源），不构成独立外部验证；真正独立的验证需 κ_L（MACE 声子）等异源标签。
- PF 在 Structure graph 上仅弱平滑（z=-3.6），在 Electronic graph 上不显著（z=-1.0）。
- 结论：**PF 最连续的空间是 Transport graph**，符合"PF 是输运性能标签"的设定。

## Eg smoothness on Structure graph

- smoothness=2.15，z=-16.7，p=0.000（显著）。
- 结构局部邻域的 band gap 显著比随机更相近（即使全局 distance correlation 仅 0.06）——结构局部邻域与电子结构存在统计连续关系。

## 说明

PF 全程使用 JARVIS-convention database-defined PF（Phase L0 已确认其 eigenvalue-pairing 敏感），
不解释为方向性 PF；PF 未参与任何图构建/调参。
