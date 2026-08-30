# PF 本征值配对歧义审计（L0-A）

## 问题
JARVIS 源码 PF_i = S_i^2 * sigma_i / 1e6，但 np.linalg.eigvals 分别作用于 S 与 sigma 张量，
二者本征值顺序未必对应同一物理方向。枚举 sigma 的 6 种 permutation，评估 PF 对配对顺序的敏感度。

## 结果

| carrier | n | ambiguity median | ambiguity mean | p90 | p95 | max | >=0.50 占比 | JARVIS rel_pos median |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| n | 806 | 0.637 | 1.005 | 1.123 | 1.509 | 223.7 | 59.6% | 0.921 |
| p | 803 | 0.473 | 1.121 | 1.066 | 1.528 | 304.7 | 47.7% | 0.917 |

歧义分档（n / p）：<0.05 占 3.3% / 3.4%；0.05–0.20 占 13.6% / 14.7%；0.20–0.50 占 23.4% / 34.2%；
>=0.50 占 59.6% / 47.7%。

## 结论
1. **PF 对 S 与 sigma 的本征值配对顺序高度敏感**：中位歧义 0.47–0.64，约 82–83% 材料歧义 > 0.20。
2. JARVIS 的 identity pairing 位于 6 种配对中**接近极大值**的位置（rel_pos≈0.92），
   这是因为 np.linalg.eigvals 对 S 与 sigma 的排序大致同向（大 |S| 配大 sigma，使 PF 偏大）。
3. 因此 PF 必须定义为 **database-defined PF under JARVIS eigenvalue pairing convention**，
   不能解释成严格 crystallographic principal-direction PF。
4. PF 仍可用于数据库内部排序 / high-PF mapping / relative comparison（所有材料同一约定），
   但**不能**过度解释方向性。这进一步支持：PF 只作 external label，不进 embedding。
