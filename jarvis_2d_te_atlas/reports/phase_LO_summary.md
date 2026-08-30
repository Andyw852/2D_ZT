# Phase L–O 执行报告：单物理视图构建、相似图建立与几何验证

## Executive Summary

### 1. Structure representation
- geometry-only SOAP（dummy species，n_max=6, l_max=6, sigma=1.0, periodic=True, r_cut=6 Å, mean-pool, L2 归一化，147 维）+ 81 元素 elemental-fraction（Hellinger）。
- 融合权重 0.5/0.5（baseline）；species-sensitive SOAP 因 81 元素维度爆炸而 SKIPPED_HIGH_DIMENSION。
- 关键修正：periodic=False → periodic=True（否则 supercell invariance 从 0.86 → 1.0）。

### 2. Graph parameters
- 所有主图统一 k=15（Structure 是使 giant=1.0 的最小 k；Electronic/Transport 沿用）。
- Structure：r_cut=6, w=0.5/0.5, distance = 0.5*d_geo_norm + 0.5*d_comp_norm。
- Electronic：RobustScaler + Euclidean；Transport：RobustScaler + Euclidean（6 个冻结 V1）。

### 3. Single-view stability
- Structure 图最稳定（giant=1.0、单分量、0 isolated；r_cut/fusion 邻域中等敏感）。
- n/p Transport 图全连通（giant=1.0、0 isolated，修复 D_sigma NaN 后）。
- Electronic 图（Eg 层与 rich 层）**天然分裂为金属/半导体两个分量**（物理正确，非缺陷），giant 0.64–0.70。

### 4. Cross-view relationship
- Structure 与 Electronic/Transport **几乎完全去相关**（距离相关 0.02–0.07；neighbor overlap ~1.5%≈随机）。
- Electronic 与 Transport 中等相关（0.32–0.38）；n-Transport 与 p-Transport 显著相关但不同（overlap 10.7% z=60，距离相关 0.42）。

### 5. PF external validation
- PF 在 n/p Transport 图上高度平滑（z≈-24～-26），在 Structure 图上仅弱平滑（z=-3.6），在 Electronic 图上不显著（z=-1.0）。
- 即 PF 最连续的物理空间是 Transport graph，验证了 transport representation 捕获了 PF 相关的基础输运物理。

---

## Final Frozen Graphs

| graph | N | k | features | distance | normalization |
|---|---:|---:|---|---|---|
| G_structure_v1 | 1103 | 15 | SOAP(r_cut=6)+elemental-fraction | 0.5*d_geo+0.5*d_comp | local-scale Gaussian, union |
| G_Eg_v1 | 1103 | 15 | Eg_optb88vdw | |Eg_i-Eg_j| | local-scale Gaussian, union |
| G_electronic_n_v1 | 678 | 15 | Eg + m_elec_median | Euclidean (RobustScaler) | local-scale Gaussian, union |
| G_electronic_p_v1 | 678 | 15 | Eg + m_hole_median | Euclidean (RobustScaler) | local-scale Gaussian, union |
| G_n_transport_v1 | 806 | 15 | 6×V1 | Euclidean (RobustScaler) | local-scale Gaussian, union |
| G_p_transport_v1 | 803 | 15 | 6×V1 | Euclidean (RobustScaler) | local-scale Gaussian, union |

sensitivity：G_electronic_joint_sensitivity、G_n/p_transport_kappa_sensitivity。

---

## 25 个问题的回答

1. vacuum axis 识别：1103/1103（100%），0 ambiguous。
2. 异常结构：无（1 个厚 slab JVASP-60475 仅影响人为 15 Å 真空测试，非数据库问题）。
3. invariance：A/B/D 全部 1.0；C 中位 1.0（49/50，1 个厚 slab 例外）。通过。
4. geometry vs composition 邻域：kNN overlap 0.12（差异很大，互补）。
5. combined structure 稳定：对 fusion weight 中等敏感（0.25/0.5=0.80, 0.5/0.75=0.75），可接受。
6. SOAP 参数：n_max=6, l_max=6, sigma=1.0, periodic=True, r_cut=6 Å, mean-pool。
7. geometry/composition weight：0.5/0.5。
8. k：15。
9. 为什么 k=15：使 Structure giant=1.0 的最小 k。
10. Structure giant fraction：1.0。
11. n-Transport 连通：giant=1.0、单分量、0 isolated（稳定）。
12. p-Transport 连通：giant=1.0、单分量、0 isolated（稳定）。
13. 影响最大的 transport feature：D_sigma（移除后 overlap 0.60）与 log_sigma_dom_geo（0.68）。
14. κ_e sensitivity vs σ 主图：距离相关 >0.96，kNN overlap ~0.72（高度相似，κ_e 可不进主模型）。
15. Eg-only vs rich Electronic：Eg 层 giant 0.64（金属/半导体二分）；rich 层 giant 0.70（同样二分），rich 层多了 m* 信息。
16. Electronic-n vs n-Transport overlap：0.020（z=1.3，不显著）。
17. Electronic-p vs p-Transport overlap：0.022（z=2.7，显著但极小）。
18. Structure vs n-Transport overlap：0.016（不显著）。
19. Structure vs p-Transport overlap：0.014（不显著）。
20. n vs p Transport 是否不同空间：是，overlap 10.7%（z=60 显著）但距离相关仅 0.42，明显不同。
21. cross-view overlap 是否都显著高于随机：否——只有 Structure-vs-Eg（z=22）、Electronic_p-vs-Transport_p（z=2.7）、n-vs-p Transport（z=60）显著；Structure-vs-Electronic/Transport 不显著。
22. PF 在 n-Transport 平滑：是（z=-24）。
23. PF 在 p-Transport 平滑：是（z=-26）。
24. PF 哪个空间最连续：Transport graph（z=-24）> Structure（z=-3.6）> Electronic（z=-1.0，不显著）。
25. 最大 Structure–Transport disagreement 材料：已输出 data/audit/high_structure_transport_disagreement_*.csv 供人工检查（本阶段仅 diagnostic，非正式 candidate）。

---

## Final Single-View Decision Table

| Layer | N | Features | Distance | k | Status | Role |
|---|---:|---|---|---:|---|---|
| Structure | 1103 | SOAP + composition | fused 0.5/0.5 | 15 | frozen | backbone |
| Eg | 1103 | Eg | |ΔEg| | 15 | frozen | electronic scalar layer |
| Electronic-n | 678 | Eg + m_elec | RobustScaler+Euclid | 15 | frozen | n electronic |
| Electronic-p | 678 | Eg + m_hole | RobustScaler+Euclid | 15 | frozen | p electronic |
| n-Transport | 806 | 6×V1 | RobustScaler+Euclid | 15 | frozen | primary |
| p-Transport | 803 | 6×V1 | RobustScaler+Euclid | 15 | frozen | primary |
| κe-n / κe-p | 806/803 | T2 | RobustScaler+Euclid | 15 | sensitivity | secondary |

---

## Ready for Phase P?

**YES** —— 所有 6 张主图 + 3 张 sensitivity 图已完成 QA / stability / connectivity / neighbor analysis / external-label diagnostics。

唯一需要在下阶段留意：Electronic View 的金属/半导体二分（2 分量），以及 Structure 与 property view 几乎去相关这一负结果——两者都是重要发现，需在 joint alignment 阶段正确处理（Structure backbone 负责桥接）。
