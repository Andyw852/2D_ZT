# 热电候选筛选状态（Step 12）

> **结论：当前数据不允许生成候选排名。** 旧版 `PF_best / κ_L_snyder` 排名已撤回。

## 前置条件验收

| 条件 | 通过 | 证据 |
|---|---:|---|
| electronic_and_kappa_same_material_id | 否 | JARVIS↔MP only reduced_formula; polymorph identity is unknown |
| same_temperature | 否 | PF at 600 K; Snyder model fixed at 300 K |
| independent_kappa_L_target | 否 | large-N target is an elastic-property formula, not measured/BTE κL |
| kappa_e_units_and_tau_validated | 否 | local JARVIS fields are not converted into a validated zT denominator |
| energy_above_hull_available | 否 | stability field is absent locally |
| unique_formula_rows | 是 | formula proxy table is de-duplicated, but this does not resolve polymorph mapping |

## 为什么撤回旧排名

旧排名把 600 K 的 JARVIS PF 与 300 K 的 Snyder 解析模型按化学式拼接，
把同一 PF 复制给多个 MP 多晶型，忽略 κe、τ、energy_above_hull，再把 PF/κL
称作 zT 排序代理。该顺序不能由 zT 公式推出，已知材料的名次也不能作为验证。

## 恢复筛选所需数据

- material_id 对齐的电子输运（PF、σ/τ、κe/τ）与载流子浓度、温度；
- 同一结构、同一温度的实验或 BTE κL；
- energy_above_hull；
- τ 敏感性和已知热电体系的盲法召回验证。