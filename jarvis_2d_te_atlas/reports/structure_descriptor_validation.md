# Structure Descriptor 验证（Phase L）

## 最终 descriptor

- Geometry-only SOAP：所有原子映射为 dummy species X，n_max=6, l_max=6, sigma=1.0, periodic=True, r_cut=6 Å, mean-pool, L2 归一化，147 维。
- Composition：81 元素 elemental-fraction vector，Hellinger distance。
- species-sensitive SOAP：SKIPPED_HIGH_DIMENSION（81 元素会得到约 48.8 万维，不可行）。

## 关键修正

- 初版用 periodic=False，导致 supercell invariance 失败（0.86）——因为 primitive cell 缺少面内周期镜像。
  改用 periodic=True（配合 pbc=[T,T,F]，dscribe 仅在面内复制、面外不复制），supercell invariance 升到 1.000000。

## Invariance tests（r_cut=6）

| test | min | median | 结论 |
|---|---:|---:|---|
| A atom permutation | 1.000000 | 1.0 | PASS |
| B slab translation | 1.000000 | 1.0 | PASS |
| C vacuum (15/20/25/30 Å) | 0.892 | 1.0 | PASS（49/50；1 个厚 slab 例外） |
| D 2×2 supercell | 1.000000 | 1.0 | PASS |

- Test C 唯一例外 JVASP-60475：slab 厚 22.88 Å，比人为的 15 Å 测试真空还厚，物理上放不下，属测试假象而非 SOAP 缺陷（其数据库真实 vacuum gap=21.74 Å）。

## 结论

SOAP 通过全部主要 invariance tests，可进入 Phase M。
