# λ 选择（P-D）

扫描 λ ∈ {0.01, 0.03, 0.1, 0.3, 1, 3, 10}，对 n/p 的 Full 与 Property 配置分别计算 anchor tension + preservation。

关键观察（n-Property，Eg+m_e+Tn）：
| λ | tension | P_transport | P_Eg | P_mass |
|---:|---:|---:|---:|---:|
| 0.1 | 0.700 | 0.536 | 0.578 | 0.623 |
| 0.3 | 0.547 | 0.607 | 0.554 | 0.622 |
| 1.0 | 0.358 | 0.555 | 0.502 | 0.507 |
| 3.0 | 0.200 | 0.368 | 0.372 | 0.266 |

- **λ=0.3 是甜点**：P_transport 最大（0.607），tension 已进入下降段，且尚未进入 λ>1 的 preservation 崩塌区。
- stable range ≈ [0.1, 1.0]（P_transport > 0.55）。
- λ>3 时 P_Eg/P_mass 崩塌（alignment 压过 preservation，即 alignment–preservation trade-off）。
- 选择 λ=0.3（transport preservation 最大 + tension 合理 + 较小 λ）。

（P_v 从未达到 0.70 理想值，最高 ~0.64；0.60 阈值下 Eg 层受金属 clique 影响 ~0.55。）
