# Joint PF 一致性检查（同源，非独立验证）

> **注意**：PF_mean = S²σ 由 nseeb×ncond 定义，而 transport 图特征本身已嵌入 S 与 σ，
> 故 PF 在联合流形上平滑是**同源一致性**，不构成对 transport 图的外部独立验证。

JARVIS-defined PF 的 log10 smoothness（1000/200 次 permutation null），在 λ=0.3 的 property 图谱 consensus kNN 图上：

| 空间 | PF smoothness z | 结论 |
|---|---:|---|
| Transport | -24.0 / -25.8 | 高度平滑 |
| **Joint** | **-6.5 / -8.4** | 显著平滑 |
| Structure | -3.6 | 弱平滑 |
| Electronic | -1.0 | 不显著 |

- PF 在 Joint 上显著平滑（z=-6.5/-8.4，p=0.000），但这是**同源一致性**（PF 与 transport 特征同源），不构成独立验证。
- 但低于纯 Transport（z=-24），因为 joint 混合了 Eg/mass 信息，稀释了 PF 信号。
