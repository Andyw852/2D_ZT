# Electronic View 特征选择（L0-H）

## 分层设计（第 51–52 节）

- **Eg electronic layer**：OptB88vdW band gap，覆盖 1103（100%）。
- **Rich electronic layer**：Eg + effective-mass spectrum，覆盖 678。
- 不因 effective mass 只覆盖 678 而删除其余 425 个材料；下一阶段 partial multi-view 决定如何组合。

## 冻结的 Electronic V1 变量

| feature | 角色 | 覆盖率 | 说明 |
|---|---|---:|---|
| Eg_optb88vdw | primary | 1103 | 0=金属 |
| m_elec_median | primary | 678 | electron 有效质量稳健中位数（面内类），0.56 m_e 中位 |
| m_hole_median | primary | 678 | hole 有效质量稳健中位数，1.14 m_e 中位 |
| m_elec_dom_geo | candidate | 678 | 两小主通道几何均值 |
| m_hole_dom_geo | candidate | 678 | 同上 |
| m_elec_spectral_ratio | candidate | 678 | 2D 质量各向异性 |
| m_hole_spectral_ratio | candidate | 678 | 同上 |
| Eg_mbj | higher-level validation | 246 | 22% 覆盖率 |
| Eg_hse | exploratory | 54 | 4.9% 覆盖率 |

## 关键决策

1. effective mass 用 **median 而非 mean**（mean 被面外类大值污染 300–575 倍）。
2. effective mass 属于 **Electronic View**（band/electronic structure state），
   不进入 Transport View（transport performance）。
3. MBJ/HSE gap 覆盖率低，不进主 Electronic View，仅作 validation/exploratory。
