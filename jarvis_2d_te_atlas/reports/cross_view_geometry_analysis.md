# Cross-View 几何分析（Phase O 核心）

## Neighbor overlap（k=15）+ random baseline

| pair | N | overlap | null mean | z | p |
|---|---:|---:|---:|---:|---:|
| Structure vs Eg | 1103 | 0.033 | 0.0125 | +22.1 | 0.000 |
| Structure vs Electronic_n | 678 | 0.016 | 0.0141 | +1.6 | 0.073 |
| Structure vs Electronic_p | 678 | 0.015 | 0.0159 | -0.6 | 0.716 |
| Structure vs Transport_n | 806 | 0.016 | 0.0141 | +1.3 | 0.098 |
| Structure vs Transport_p | 803 | 0.014 | 0.0140 | -0.1 | 0.558 |
| Electronic_n vs Transport_n | 675 | 0.020 | 0.0183 | +1.3 | 0.102 |
| Electronic_p vs Transport_p | 674 | 0.022 | 0.0185 | +2.7 | 0.005 |
| Transport_n vs Transport_p | 802 | 0.107 | 0.0193 | +60.6 | 0.000 |

## Distance correlation（Spearman）

| pair | rho |
|---|---:|
| Structure vs Eg | 0.056 |
| Structure vs Electronic_n/p | 0.049 / 0.066 |
| Structure vs Transport_n/p | 0.026 / 0.024 |
| Electronic_n/p vs Transport_n/p | 0.380 / 0.320 |
| Transport_n vs Transport_p | 0.417 |
| Electronic_n vs Electronic_p | 0.743 |

## 核心结论

1. **Structure 与所有 property view 几乎完全去相关**（距离相关 0.02–0.07，neighbor overlap ~1.5% 与随机相当）。
   即：结构相似 **不能** 推出电子/输运相似。这是本研究最重要的负结果之一。
2. **n-Transport 与 p-Transport 是显著相关但明显不同的两个空间**（overlap 10.7% >> 随机 1.9%，z=60；距离相关 0.42），支持 n/p 分开建模。
3. **Electronic 与 Transport 有中等距离相关**（0.32–0.38），符合物理预期（Eg+m* 部分决定输运）。
4. 并非所有 cross-view overlap 都显著高于随机：Structure-vs-property 大多不显著（p=0.07–0.72）。
