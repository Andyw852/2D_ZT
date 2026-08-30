# 全体材料总体数据图：结果与边界

## 核心修正

原模型图只代表一个基准参数组合，不能解释为材料总体分布。本目录的主图现已改为：每个关系使用它所需的最大成对完整样本集；真实 ZT 面板每个实验样品只出现一次，其他物性在同一样品的峰值 ZT 温度配对。

## 数据覆盖

| 性质 | 直接实验样品数 |
|---|---:|
| Seebeck | 33,950 |
| total thermal conductivity | 22,259 |
| ZT | 20,037 |
| electrical conductivity | 17,059 |
| power factor | 16,807 |
| lattice thermal conductivity | 6,371 |
| carrier concentration | 1,763 |
| mobility | 1,524 |
| electronic thermal conductivity | 737 |

峰值 ZT 配对数：

| 关系 | 同样品配对数 |
|---|---:|
| ZT + n | 440 |
| ZT + S | 12,070 |
| ZT + sigma | 11,738 |
| ZT + PF | 11,937 |
| ZT + kL | 3,649 |
| ZT + kTotal | 11,859 |

## 真实 ZT 总体关系

| 横轴 | N | Spearman rho(ZT,x) |
|---|---:|---:|
| n | 440 | -0.178 |
| |S| | 12,070 | +0.364 |
| sigma | 11,738 | +0.282 |
| kappaL | 3,649 | -0.438 |

散点密度和分箱中位数才是主要读图对象，不能用单条 SPB 曲线代替总体数据。由于样品来自不同化学家族、温度、制备状态和论文，相关性是描述性的，不是受控因果效应。
电导配对优先使用直接电导曲线；样品只有电阻率曲线时使用 sigma=1/rho 补充，因此电导配对数可以高于直接电导覆盖数。

## 结构数据能回答什么

- 相对密度可解析样品：2,201；可近似得到孔隙率，但许多记录是阈值或区间。
- 晶粒尺寸可解析样品：638；用于分档，不把文本中值当高精度测量。
- 弹性—实验 κL 有效跨库配对：79；可检验 sqrt(B/rho) 声速代理，但存在化学式多晶型歧义。
- JARVIS 2D 几何—输运记录：1,609 行；用于形状代理与 S/sigma，不是真实 ZT。
- 当前数据没有显式褶皱幅度/波长、孔洞形貌、声子群速度或完整声学/光学支标签，因此不能把这些变量画成经验总体关系；模型图只作为待验证假设保留。

## 图清单

1. `01_data_coverage.png`：每种性质及每个配对关系到底用了多少样品；
2. `02_experimental_ZT_global.png`：真实 ZT 与 n、S、sigma、kappaL 的全体可配对样品；
3. `03_electronic_all_available.png`：不要求 ZT 的 Pisarenko、电导、迁移率、PF 总体数据；
4. `04_thermal_all_available.png`：kappaL 温度带、总/晶格热导、电子热导及声速代理；
5. `05_structure_metadata_experimental.png`：真实 ZT 与孔隙率、晶粒、样品形态；
6. `06_jarvis2d_shape_transport.png`：二维原子几何代理与 S/sigma。

每个面板的来源、样本数和温度规则见 `panel_manifest.csv`。
