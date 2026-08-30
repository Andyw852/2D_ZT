# Phase P–S 执行报告：不完备多视图对齐、公共潜在空间验证与统一 Transport Atlas

## Executive Summary

### A. n-Atlas
**SUPPORTED**（Eg + electron mass + n-Transport 属性图谱，moderate preservation）

### B. p-Atlas
**SUPPORTED**（Eg + hole mass + p-Transport 属性图谱，moderate preservation）

### C. Structure
**AUXILIARY ONLY**（不进公共空间）

### D. λ
- lambda_n = **0.3**，lambda_p = **0.3**
- stable range ≈ **[0.1, 1.0]**（P_transport 在该区间 > 0.55，tension 已进入下降段，λ>1 后 preservation 崩塌）

### E. Main Scientific Finding
> **Electronic–Transport 确实形成比 Structure–Transport 明显更强的共同材料几何。**
- Electronic vs Transport distance correlation = 0.32–0.38（中等），joint preservation P_transport = 0.61–0.64。
- Structure vs Transport distance correlation = 0.02–0.07（≈随机），joint preservation P_structure = 0.13–0.51（差）。
- 结论：**Structure 与 property views 是本质不同的材料空间**；二维材料的结构相似性不能预测其电子/输运相似性。公共材料表示只对 property views（Electronic + Transport）成立，Structure 只能作为 linked auxiliary view。

---

## Final Architecture

\`\`\`text
FINAL N-TYPE ATLAS V1                    FINAL P-TYPE ATLAS V1
Consensus Identity Node C_i              Consensus Identity Node C_i
        |                                       |
        +--- Eg (Eg_optb88vdw)                  +--- Eg
        +--- electron effective mass            +--- hole effective mass
        +--- n-Transport (6×V1)                 +--- p-Transport (6×V1)

Structure: linked auxiliary view        Structure: linked auxiliary view
PF: external label                      PF: external label
kappa_e: sensitivity only               kappa_e: sensitivity only
\`\`\`

## Frozen Joint Models

| model | views | N | lambda | alpha | latent dim | P_transport | P_Eg | P_mass | random-anchor |
|---|---|---:|---|---:|---:|---:|---:|---|
| n-Property | Eg + m_e + Tn | 1103 (core 806) | 0.3 | 1 | 20 | 0.607 | 0.554 | 0.622 | z=4.4 (passed) |
| p-Property | Eg + m_h + Tp | 1103 (core 803) | 0.3 | 1 | 20 | 0.637 | 0.558 | 0.603 | z=5.8 (passed) |

（n-Full/p-Full 含 Structure，P_structure 仅 0.13–0.51，且 Structure 与 property 去相关 → 不采用。）

## PF External Validation（JARVIS-defined PF, log10 平滑度 z-score）

| 空间 | z |
|---|---:|
| Structure | -3.6 |
| Electronic | -1.0（不显著） |
| Transport | -24.0 / -25.8 |
| **Joint** | **-6.5 / -8.4** |

- PF 在 Joint 上显著平滑（z=-6.5/-8.4），但低于纯 Transport（z=-24）。
- 即 Joint 空间捕获了部分 PF 结构（多于 Structure/Electronic，少于纯 Transport），符合"PF 是输运性能标签"。

## Cross-View Physics

- Structure ↔ Electronic / Transport：distance correlation 0.02–0.07，neighbor overlap ~1.5%（≈随机）→ 几乎独立。
- Electronic ↔ Transport：0.32–0.38（中等）→ 共享部分几何。
- n ↔ p：consensus distance Spearman 0.766，kNN(15) overlap 0.495 → 相关但明显不同，需分别建 Atlas。

## Stable Anomalous Materials（cross-view transport anomaly）

- n 高 T_transport：JVASP-14456, JVASP-27864, JVASP-7033, JVASP-27853, JVASP-20012。
- p 高 T_transport：JVASP-13619, JVASP-5947, JVASP-6601, JVASP-6805, JVASP-786。
- （仅称 anomalous materials，不称 superlattice candidates。完整 top-50 见 data/processed/joint_tension_*.csv。）

## Final Status Variables

\`\`\`text
N_ATLAS_SUPPORTED                = True
P_ATLAS_SUPPORTED                = True
STRUCTURE_IN_N_COMMON_SPACE      = False
STRUCTURE_IN_P_COMMON_SPACE      = False
COMBINED_NP_ATLAS_SUPPORTED      = False (n/p 明显不同，分别建 Atlas)
UNIFIED_MANIFOLD_SUPPORTED       = True (property-only；Structure 为 linked)
RANDOM_ANCHOR_CONTROL_PASSED     = True (preservation z=4.4/5.8)
DUPLICATE_SENSITIVITY_PASSED     = True (12 组精确重复；真正重跑 collapsed 流形，kNN overlap n=0.96/p=0.92)
\`\`\`

## Ready for Superlattice Parent Discovery?

**YES** —— 但注意：由于 Structure 不在公共空间，下一阶段应依据 linked-view 的 structure/electronic/transport 距离 + cross-view disagreement 筛选材料对（见 phase 46 的 Linked 分支），而非单一 joint 空间坐标。

---

## 28 个问题的答案

1. n joint manifold：SUPPORTED（moderate，P_transport=0.61）。
2. p joint manifold：SUPPORTED（moderate，P_transport=0.64）。
3. lambda_n = 0.3。
4. lambda_p = 0.3。
5. stable range ≈ [0.1, 1.0]。
6. Structure 不进 n-atlas（auxiliary）。
7. Structure 不进 p-atlas（auxiliary）。
8. n-Property 优于 n-Full（结构引入不提升且 P_structure 低）。
9. p-Property 优于 p-Full。
10. Transport preservation：n 0.607 / p 0.637。
11. Eg preservation：n 0.554 / p 0.558。
12. Mass preservation：n 0.622 / p 0.603。
13. Structure preservation：0.13–0.51（差）。
14. real vs random anchors：preservation z=4.4（n）/ 5.8（p），显著优于随机。
15. 影响最大的 layer：Eg（drop 后 consensus overlap 0.217），其次 transport（0.392），mass（0.575）。
16. duplicate collapse：12 组，影响小，结论不变。
17. n/p Atlas 属于不同公共 geometry：是（consensus distance 相关 0.766，kNN overlap 0.495）。
18. Electronic–Transport 共享几何强于 Structure–Transport：是（0.32–0.38 vs 0.02–0.07）。
19. metal/semiconductor boundary：仍存在（Eg 是最强 axis，397 metal / 706 semi）。
20. PF 在 Joint 显著平滑：是（z=-6.5/-8.4）。
21. Joint > Structure：是（-6.5 vs -3.6）。
22. Joint > Electronic：是（-6.5 vs -1.0）。
23. Joint > Transport：否（-6.5 < -24，Joint 稀释了 PF 信号）。
24. T_structure 最大材料：见 joint_tension 表（Structure 未入公共空间，T_structure 反映其与 property 表示不一致）。
25. T_transport 最大材料：n=JVASP-14456 等，p=JVASP-13619 等。
26. high-tension 在参数变化下稳定性：需 λ 重扫确认（本阶段未全量重扫，标注待验证）。
27. Combined n+p Atlas：不保留（n/p 明显不同）。
28. Unified 还是 Linked：**Linked Multi-View Atlas（property 部分统一，Structure linked）**。
