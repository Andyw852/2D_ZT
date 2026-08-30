# 2D Thermoelectric Materials Atlas (二维热电材料多物理统一图谱)

利用 JARVIS curated dft_2d 数据库，在不进行新第一性原理计算、不预测/填补缺失性质的前提下，
分别建立结构、电子结构、热电输运性质对应的材料相似性空间，并通过不完备多视图流形对齐
（JID 身份锚点 + 多层图谱嵌入）映射到统一的二维材料公共潜在空间。

## 当前阶段（第一轮）：Phase A–F

- Phase A: 建立 Python 环境 (te_manifold, Python 3.11)
- Phase B: 验证 dft_2d 真实官方下载源
- Phase C: 下载 dft_2d
- Phase D: 保存原始快照 + SHA256
- Phase E: Schema Audit（全部字段）
- Phase F: Property Coverage Audit
- STOP —— 得到真实 Property Coverage Matrix 之前不进行任何流形建模

## 核心约束

- 仅使用 JARVIS curated dft_2d 一个数据库
- 缺失值绝不人为设置为 0 或插值/预测
- 有该性质才建立该性质节点/层
- n 型 / p 型 Seebeck 与 PF 禁止平均，分别研究
- ZT 只作为外部验证标签，不进入流形构建特征（target leakage 防护）

## 现有数据的发现能力验证

`scripts/44_masked_transport_manifold_validation.py` 对组成/结构、电子结构和 n/p 电输运执行
隐藏视图验证：每一折先删除测试材料的全部输运特征，再用训练材料的输运图帮助流形对齐，
以数据库定义的 PF 作独立测试标签。结果见
`reports/masked_transport_manifold_validation.md`。该验证不使用声子、晶格热导率、稳定性或
可合成性数据，结论限定为电输运高-PF 候选富集，不解释为真实 ZT 预测。

## 目录结构

见 scripts/ 下的 00–16 脚本。
