# Structure Graph 选择（Phase M）

## 参数选择

- SOAP r_cut=6 Å（r_cut 4-6 overlap 0.56、6-8 overlap 0.70，6 是合理中间值）。
- fusion weight = 0.5/0.5（预注册 baseline；0.25/0.5 overlap 0.80、0.5/0.75 0.75，对权重中等敏感）。
- k=15（使 giant component 达到 1.0、单连通分量的最小 k）。

## k 扫描（r_cut=6, w=0.5）

| k | giant | components | isolated |
|---:|---:|---:|---:|
| 10 | 0.990 | 2 | 0 |
| 15 | 1.000 | 1 | 0 |
| 20 | 1.000 | 1 | 0 |

## 结构邻域由什么驱动？

- geometry-only vs composition-only 的 kNN(20) overlap = 0.12（二者给出非常不同的邻域）。
- combined(0.5/0.5) 与 composition-only overlap 0.55、与 geometry-only 0.27 → combined 图更偏 composition 驱动。

## Near-duplicates

- d_struct 最低的前若干对都是 d=0.0 的同 formula 精确重复（不同 JID），如 NbS2(JVASP-5989 vs 5992)、MoSe2、Bi2Pt、FeSe、BN 等。
- 这是 dft_2d 库中同结构不同 JID 的多晶型/重复条目，后续多视图对齐需注意。

## 冻结

G_structure_v1：N=1103，k=15，giant=1.0，单分量，0 isolated，mean degree 21.3。
