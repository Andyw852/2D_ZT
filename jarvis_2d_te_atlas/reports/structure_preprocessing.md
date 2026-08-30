# 二维结构预处理（Phase L）

## vacuum axis 识别

- 全部 1103 个结构都成功识别出 vacuum axis（0 个 ambiguous，confidence min=3.8 >> 1.5）。
- vacuum axis 分布：1102 个为第 3 轴（z），1 个为第 1 轴。
- vacuum_confidence median=12.0（即最大空腔 gap 约是次大 gap 的 12 倍），识别非常明确。

## 标准化

- 将 vacuum axis 重排为第 3 方向，沿该方向把 slab 质心居中，pbc=[True, True, False]。
- 保存到 data/processed/standardized_2d_structures.parquet（1103 条），原始数据未覆盖。

## 结论

无异常二维结构；无需 STOP。SOAP 可以继续。
