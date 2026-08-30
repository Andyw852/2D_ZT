# Seebeck 符号一致性审计（L0-C）

## 符号一致性分类

| class | n-type 含义 | n 数量 | p-type 含义 | p 数量 |
|---|---|---:|---|---:|
| 3/3 | N3 全负 | 564 (70.0%) | P3 全正 | 601 (74.8%) |
| 2/3 | N2 | 131 (16.3%) | P2 | 86 (10.7%) |
| 1/3 | N1 | 54 (6.7%) | P1 | 54 (6.7%) |
| 0/3 | N0 | 57 (7.1%) | P0 | 62 (7.7%) |

sign_fraction 均值：n = 0.83，p = 0.84（平均 83–84% 的本征值符合预期符号）。

## mean vs median 稳健性

- sign(S_mean) == majority_sign 的比例：n 92.8%，p 93.2%。
- sign(S_median) == majority_sign 的比例：**n 100%，p 100%**。
- 结论：**S_median 比 S_mean 更稳健**，median 天然不受单个符号异常本征值影响，始终与多数本征值符号一致。

## 符号异常与带隙的关系（假设检验 H1/H2/H3）

| class | n gap median (eV) | n gap mean | p gap median (eV) |
|---|---:|---:|---:|
| 3/3 (正常) | 1.447 (N3) | 1.641 | 1.517 (P3) |
| 2/3 | 0.401 (N2) | 1.002 | 0.000 (P2) |
| 1/3 | 0.000 (N1) | 0.000 | 0.000 (P1) |
| 0/3 (严重违例) | 0.000 (N0) | 0.002 | 0.000 (P0) |

- bad-sign（N0/N1 或 P0/P1）gap median = 0.000 eV（金属），good-sign（N2/N3 或 P2/P3）gap median = 1.34–1.36 eV。
- Mann-Whitney U 检验 p = 3.1e-47（n）/ 4.0e-47（p），**极显著**。

## 结论（措辞：associated with，非 caused by）

1. 约 15% 的 Seebeck 符号违例**高度集中在 zero/small band gap（金属）区域**（p~1e-47）。
   这支持 **H2（small-gap / bipolar transport）**，而非单纯的 H1（数值不稳定）。
2. 符号违例材料未删除，JID 保留在 data/audit/seebeck_sign_violations.csv 供人工复核。
3. 处理原则：n/p Transport 主模型使用 **S_median（稳健）** 而非 S_mean；
   并额外保留 S_sign_fraction 作为 bipolar/metallicity 的物理标志（与 Eg 强相关，但不等于 Eg）。
