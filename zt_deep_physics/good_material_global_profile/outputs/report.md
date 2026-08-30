# 高性能热电材料相对全局的必要特征

## 结论先行

本数据集共有 **20,037** 个具有有效峰值 ZT 的实验样品。主分析把全局前 10% 定义为高性能组：`ZT_peak >= 1.227`，共 **2,004** 个样品；另有 **4,002** 个样品满足 `ZT >= 1`。

数据支持的核心不是一个孤立阈值，而是同时满足：**中等到较高的 |S|、足够高的电导/功率因子，以及较低的晶格或总热导**。载流子浓度给出可用工作窗，但单独区分高 ZT 的能力有限；迁移率和结构元数据当前覆盖不足或分离较弱，不能设成普适硬门槛。

下表的 `P10-P90` 覆盖 80% 高性能配对样品，适合作为软筛选盒；它不是数学意义上的必要条件。

| 性质 | 高性能组 N / 覆盖率 | 全局中位数 | 高性能典型 P10-P90 | 高性能中位数 | 相对全局 | 证据 |
|---|---:|---:|---:|---:|---|---|
| Peak-ZT temperature | 2,004 / 100.0% | 723 K | 401-900 K | 768 K | 中等偏高 (effect=0.12) | strong |
| Carrier concentration | 84 / 4.2% | 3.90e+19 cm^-3 | 5.79e+18-6.14e+20 cm^-3 | 3.25e+19 cm^-3 | 中等偏低 (effect=0.13) | limited |
| Absolute Seebeck | 1,231 / 61.4% | 189 uV/K | 186-299 uV/K | 230 uV/K | 明显偏高 (effect=0.45) | strong |
| Electrical conductivity | 1,217 / 60.7% | 41342 S/m | 15862-1.15e+05 S/m | 52809 S/m | 中等偏高 (effect=0.13) | strong |
| Power factor | 1,261 / 62.9% | 1.26 mW/mK^2 | 0.998-4.69 mW/mK^2 | 2.75 mW/mK^2 | 明显偏高 (effect=0.51) | strong |
| Lattice thermal conductivity | 550 / 27.4% | 0.814 W/mK | 0.298-0.999 W/mK | 0.591 W/mK | 明显偏低 (effect=0.43) | strong |
| Total thermal conductivity | 1,238 / 61.8% | 1.62 W/mK | 0.578-2.52 W/mK | 1.10 W/mK | 明显偏低 (effect=0.32) | strong |
| Electronic thermal conductivity | 36 / 1.8% | 0.723 W/mK | 0.209-1.13 W/mK | 0.455 W/mK | 中等偏低 (effect=0.18) | limited |
| Mobility | 92 / 4.6% | 39.6 cm^2/Vs | 2.24-139 cm^2/Vs | 32.9 cm^2/Vs | 中等偏低 (effect=0.20) | limited |
| Relative density | 66 / 3.3% | 95.0 % | 95.0-98.0 % | 97.0 % | 明显偏高 (effect=0.59) | limited |
| Porosity fraction | 66 / 3.3% | 0.050 fraction | 0.020-0.050 fraction | 0.030 fraction | 明显偏低 (effect=0.59) | limited |
| Grain size | 22 / 1.1% | 1.00 um | 0.033-1.00 um | 1.00 um | 中等偏低 (effect=0.29) | limited |

## 哪些特征最接近经验上的必要条件

- **Power factor**：高性能组典型范围 `0.998-4.69 mW/mK^2`；明显偏高，高性能样品覆盖率 62.9%。
- **Absolute Seebeck**：高性能组典型范围 `186-299 uV/K`；明显偏高，高性能样品覆盖率 61.4%。
- **Lattice thermal conductivity**：高性能组典型范围 `0.298-0.999 W/mK`；明显偏低，高性能样品覆盖率 27.4%。
- **Total thermal conductivity**：高性能组典型范围 `0.578-2.52 W/mK`；明显偏低，高性能样品覆盖率 61.8%。
- **Electrical conductivity**：高性能组典型范围 `15862-1.15e+05 S/m`；中等偏高，高性能样品覆盖率 60.7%。
- **Peak-ZT temperature**：高性能组典型范围 `401-900 K`；中等偏高，高性能样品覆盖率 100.0%。

### 可直接用于初筛的单侧软约束

下列界限各自保留该性质有数据的高性能样品约 90%；它们比双侧典型窗口更接近“必要条件”的操作定义。

| 性质 | 软约束 | 配对高性能覆盖率 | 温度分层后效应 |
|---|---:|---:|---:|
| Power factor | `>= 0.998 mW/mK^2` | 62.9% | 0.52 |
| Absolute Seebeck | `>= 186 uV/K` | 61.4% | 0.45 |
| Lattice thermal conductivity | `<= 0.999 W/mK` | 27.4% | 0.45 |
| Total thermal conductivity | `<= 2.52 W/mK` | 61.8% | 0.36 |
| Electrical conductivity | `>= 15862 S/m` | 60.7% | 0.12 |

### 联合规则的保留率与富集

联合规则只在所需性质同时有数据的样品中评估；`保留率`越高越接近必要，`通过后高性能率/富集`越高越适合筛选。

| 规则 | 可评估 N | 全局通过率 | 高性能保留率 | 通过后高性能率 | 富集 |
|---|---:|---:|---:|---:|---:|
| PF soft floor | 11,928 | 57.7% | 90.0% | 16.5% | 1.56x |
| kappa_total soft ceiling | 11,859 | 66.2% | 90.0% | 14.2% | 1.36x |
| kappaL soft ceiling | 3,649 | 60.5% | 90.0% | 22.4% | 1.49x |
| PF floor + kappa_total ceiling | 9,093 | 38.7% | 82.8% | 23.4% | 2.14x |
| S-sigma window + kappa_total ceiling | 8,040 | 23.7% | 70.0% | 32.8% | 2.96x |

## 富集而不是因果

每个性质按全局十分位分箱后，最高富集区如下。富集倍数以该性质的成对完整样本为基线，因此不会把缺失率差异误写成整体 10% 基线。

| 性质 | 最富集全局分位箱 | 数值范围 | 高性能率 | 相对基线富集 |
|---|---:|---:|---:|---:|
| Peak-ZT temperature | D6 | 723-772 K | 16.8% | 1.68x |
| Carrier concentration | D5 | 2.75e+19-3.88e+19 cm^-3 | 43.2% | 2.26x |
| Absolute Seebeck | D9 | 253-299 uV/K | 24.2% | 2.37x |
| Electrical conductivity | D7 | 57303-74705 S/m | 16.6% | 1.60x |
| Power factor | D10 | 3.74-188 mW/mK^2 | 26.8% | 2.54x |
| Lattice thermal conductivity | D3 | 0.495-0.599 W/mK | 31.0% | 2.05x |
| Total thermal conductivity | D3 | 0.835-1.04 W/mK | 19.2% | 1.84x |
| Electronic thermal conductivity | D3 | 0.233-0.353 W/mK | 23.3% | 2.78x |
| Mobility | D2 | 2.75-7.63 cm^2/Vs | 42.5% | 1.85x |
| Relative density | D6 | 96.5-100 % | 29.9% | 5.33x |
| Porosity fraction | D1 | 0.000-0.035 fraction | 26.9% | 4.80x |
| Grain size | D1 | 0.013-0.053 um | 20.0% | 1.78x |

## 结构与类别变量

相对密度、孔隙率和晶粒尺寸来自文本元数据解析，覆盖远低于电子输运数据；样品形态和材料家族的富集还会混合研究热点、温区和发表选择，故只用于提出假设。

| 类别 | 取值 | N | 高性能占比 | 富集 |
|---|---|---:|---:|---:|
| Material family | SnSe | 86 | 52.3% | 5.09x |
| Material family | Telluride | 601 | 33.9% | 3.30x |
| Sample form | Polycrystal | 449 | 33.0% | 3.20x |
| Material family | PbTe | 641 | 28.7% | 2.79x |
| Sample form | Pellet/compact | 237 | 24.9% | 2.42x |
| Material family | Alloy | 214 | 15.0% | 1.45x |
| Sample form | Powder | 48 | 14.6% | 1.42x |
| Material family | Antimonide | 478 | 13.2% | 1.28x |

## 与项目深层物理模型的对照

实验全局表没有完整的有效质量、谷简并、形变势、二维弹性常数、群速度和声子寿命。下表来自项目原有 400 个透明参数情景中、各自优化掺杂后 ZT 前 10% 的 `P10-P90`，只能作为计算先验，不能称为相对全局的经验必要范围。

| 深层变量 | 模型情景 P10-P90 | 中位数 | 证据性质 |
|---|---:|---:|---|
| DOS effective mass | 0.473-1.67 m_e | 1.11 m_e | 情景模型先验 |
| Conductivity effective mass | 0.144-0.505 m_e | 0.239 m_e | 情景模型先验 |
| Valley degeneracy | 2.90-6.00 count | 5.00 count | 情景模型先验 |
| 2D elastic modulus | 41.2-104 N/m | 78.1 N/m | 情景模型先验 |
| Deformation potential | 3.28-6.25 eV | 3.95 eV | 情景模型先验 |
| Porosity | 0.022-0.258 fraction | 0.136 fraction | 情景模型先验 |
| Phonon group velocity | 1109-2540 m/s | 1766 m/s | 情景模型先验 |
| Effective phonon lifetime | 0.082-0.393 ps | 0.203 ps | 情景模型先验 |
| Model lattice thermal conductivity | 0.114-0.755 W/mK | 0.321 W/mK | 情景模型先验 |

特别需要保留的模型—实验张力：模型高表现情景允许较宽的孔隙率窗口，但实验结构元数据中的高 ZT 配对样品主要集中在高相对密度、低孔隙率。这说明‘用孔隙降 kappaL’不是免费的收益；真实材料中的电导损失、连通性和样品制备必须同时验证。

模型与实验对 kappaL 的有利区间有明显重叠，但实验窗口更宽。对缺失的深层电子/声子变量，应把模型范围用于安排 DFT/声子计算优先级，而非直接淘汰材料。

## 如何使用这些范围

1. 初筛时先用 `screening_ranges.csv` 的 `P10-P90` 作为软筛选窗口，不要把边界外样品直接判死刑；
2. 优先联合判断 `PF` 与 `kappa_total`（或 `kappaL + kappae`），因为 ZT 本身是这些量的组合；
3. 用 `P05-P95` 检查候选是否明显偏离已知高性能经验域；
4. 二维材料还需重新核对有效厚度定义，体单位的 sigma、kappa 和载流子浓度会随厚度约定改变；
5. 对最终候选补充带隙、有效质量、谷简并、形变势、声子稳定性和三阶力常数。当前实验表并不完整包含这些深层变量。

## 证据边界

- 这是跨论文、跨材料家族、跨温度和跨制备状态的观察性比较，只能给出经验约束与优先级，不能证明单变量因果；
- 峰值 ZT 温度本身是测试条件，不是材料固有常数；高温富集也受到测量温区覆盖影响；
- 数值表同时给出温度五分位分层和材料家族内的秩效应；它们用于检查方向稳健性，不等同于完整因果校正；
- 每个性质使用各自的最大配对样本，表中的 N 不同；严禁把小样本结构范围与万级电子输运范围视为同等证据；
- PF 优先采用直接曲线，缺失时以同温度 `S^2 sigma` 补充；相对密度和晶粒尺寸是文本解析值；
- P10-P90 按定义仍允许 20% 已知高性能样品位于区间外，因此这里称为软范围而非硬必要条件。

## 可复现性

输入：`/home/wangchao/work_wc/2D_ZT/zt_deep_physics/empirical/outputs/experimental_ZT_with_structure_metadata.csv`

运行 `analyze_good_materials.py` 可重建全部表格、图和本报告。
