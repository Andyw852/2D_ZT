# κe–κL 双通道与多视图分析（审计修正版）

本目录是对原 `kappaL-multiview-verify` 管线的可复现重构与二次审计。当前版本修复了
共同样本 kNN、ID 并列破平局、实验多晶型复制、伪条件互信息、偏 Spearman、标定曲线倒推、
SOAP 默认标准化和不满足前置条件仍输出候选排名等问题。

## 正确结论（2026-08-28）

### 1. κe 的大样本“空间关系”主要是同一输运计算内部的关系

JARVIS 的电子输运字段来自 600 K、固定掺杂浓度的 CRTA/BoltzTraP 计算，并按规范化化学式
与 MP 汇总表连接。对 7,099–7,104 个唯一化学式：

| 性质 → κe | n 型 Spearman ρ | p 型 Spearman ρ | 解释 |
|---|---:|---:|---|
| 电导率 σ | +0.987 | +0.994 | 同一输运计算，Wiedemann–Franz/CRTA 方向；不是独立发现 |
| 带隙 Eg | −0.772 | −0.779 | 金属/小带隙体系在该固定条件下占主导 |
| |S| | −0.770 | −0.779 | 同一输运计算中的权衡 |
| 密度 | +0.558 | +0.564 | 跨库化学式层面的组成关联，不是因果 |

在目标未参与构造的性质空间里，加入 σ 后，空间距离与 κe 距离的 Spearman 从
0.58→0.70（n）和 0.63→0.76（p）；10-NN 富集从 3.26→16.20 和 2.89→18.28。
因此 κe 的局部邻域几乎由 σ 锚定，这与定义和计算流程一致。

### 2. 大样本 Snyder κL 的强关系主要是模型内生关系

MP 的 `snyder_acoustic` 在这里统一标作 **300 K Snyder 解析 κL 模型**，不是实验 κL，也不是
AFLOW-AGL/BTE 真值。它与剪切模量、体模量、Debye 温度的 ρ 分别为 +0.974、+0.860、
+0.852；这些变量参与或高度接近解析模型输入，不能作为新物理发现。晶格描述符空间与该模型
的距离相关为 0.765，10-NN 富集 8.72，同样主要是正对照。

### 3. 实验 κL 确实保留了弹性/声子性质信号，但样本仍小

Starrydata2 先在化学式层面聚合，不再把一个实验值复制到所有 MP 多晶型。与 MP 描述符有
96 个化学式交集；其中 70 个同时有 JARVIS 电子输运：

| 性质 → 实验 κL | N | Spearman ρ | formula-bootstrap 95% CI |
|---|---:|---:|---:|
| Debye 温度 | 96 | +0.566 | [+0.391, +0.710] |
| 剪切模量 | 96 | +0.533 | [+0.373, +0.671] |
| 体模量 | 96 | +0.427 | [+0.248, +0.595] |
| 密度 | 96 | −0.116 | [−0.321, +0.082] |

只保留“MP 中单一多晶型”的 53 个化学式时，前三项仍为 +0.498、+0.496、+0.460。
这是当前最可信的物理方向性结果，但实验温度仍混合、化学式连接不能识别具体实验晶相。

### 4. 当前数据不能证明 κe 与真实 κL 解耦

| κe ↔ κL | N | n 型 ρ [95% CI] | p 型 ρ [95% CI] | 判定 |
|---|---:|---:|---:|---|
| Snyder 模型代理 | 7,104 | +0.261 [+0.236,+0.284] | +0.266 [+0.242,+0.291] | 有弱相关，但跨温度、跨库且目标是模型 |
| 实验 κL | 70 | +0.049 [−0.216,+0.308] | +0.113 [−0.113,+0.368] | 未检出；区间过宽，不能等同于独立/正交 |

实验队列的电子性质空间与 κL 的全局距离相关约 −0.04～+0.02；10-NN 富集约 1.17～1.29。
晶格描述符空间则为距离相关 +0.235、10-NN 富集 1.47。后者方向正确但 N=70，局部效应弱。

旧版“PF⊥κL、重要性夹角 60°，所以部分解耦”的结论已撤回。当前唯一可保留的是：公式层面
PF–Snyder κL 的 ρ 约 −0.046（n）/−0.044（p），它只是 600 K JARVIS 与 300 K MP 模型的
探索性拼接，不能回答真实双通道独立调控。

### 5. Q1/Q2/Q3 的最终状态

- **Q1 结构增量信息：对真实 κL 尚无定论。** Snyder/Clarke 模型上，给已有弹性块加入 SOAP
  的配对 ΔR² 分别为 −0.00091 和 +0.00056，主要说明模型可由自身输入重构。实验 κL 的
  ΔR²=+0.057±0.067，五折范围 −0.052～+0.128，跨零且 N=59。
- **Q2 双通道：未辨识。** 实验 κe–κL 全样本相关跨零；最严格“单 MP 多晶型且单实验曲线”
  子集只有 N=22，结果对筛选敏感，仍不能证明耦合或解耦。
- **Q3 候选：阻断。** 缺少同 material_id、同温度的 PF/κe/κL、有效 τ/单位标定和
  `energy_above_hull`，因此 `candidates.csv` 保持为空，不输出伪 zT 排名。

### 6. 旧“双通道簇”图已降级，主框架改为不对称预测

`structure_space_channel_clusters.png` 等旧图不再支持“两个好材料簇”或“共享的有利结构带”。
原因不是图画得不够突出，而是标签和推断路径不成立：

- 低 Snyder κL 由平均质量、原子体积、Debye 温度、Grüneisen 参数等解析输入生成；纯几何/
  组分空间中的高可分性主要是模型内生正对照；
- JARVIS PF、κe 是 600 K、固定载流子浓度的 CRTA 输出，缺少材料相关 τ；它们不能与 300 K
  Snyder κL 合成绝对 zT 或品质因子；
- 前 5% × 后 5% 在 N≈6,553 时随机交集期望仅约 16 个，二维投影中的少量高亮点不足以定义簇；
- “唯一化学式”筛选排除了多晶型，存在系统性选择偏差。

主问题现改为：在按 chemical system 分组的相同 CV 折上，`composition + geometry` 相对
`composition-only` 能否稳定提高 `Eg`、`m*`、`ε` 的 OOF 预测。残差解释为表示/模型/协议/
匹配缺口，不再称为 `electronic-private`。只有几何增量通过预注册门槛后，才构建电子感知的结构
表示；最终热电图必须使用同温度的 `log μW` 与独立 `log κL`，叠加 `log B` 等值线和 Pareto
前沿。完整方案见 `reports/03_asymmetric_structure_to_transport_plan.md`。

### 7. 判决性实验：可以继续统一表示，但不能继续候选排名

现已直接在 55,723 条 JARVIS 原生结构上完成 chemical-system 分组 CV，多晶型不再因
“唯一化学式”筛选被删除。composition-blind SOAP 相对 162 维组成基线的增量为：

| 目标 | Ridge ΔR² | ExtraTrees ΔR² | 当前判定 |
|---|---:|---:|---|
| 正带隙大小 | +0.019 | +0.021 | 两模型 CI 均不跨 0 |
| 电子有效质量 | −0.022 | +0.021 | 模型依赖 |
| 空穴有效质量 | −0.006 | +0.014 | 模型依赖 |
| 介电响应 | +0.037 | +0.029 | 两模型 CI 均不跨 0 |

这满足“至少两个电子目标在两类模型上稳定”的进入门槛，说明几何在组分之外确有小而可复现的
电子信息；但有效质量映射仍弱。Snyder 模型输入代理的 OOF R²=0.994–0.999，确认旧低 κL
可分性主要是泄漏正对照；实验 κL 只有 N=59，几何 ΔR² 的 CI 仍跨 0。

结果图与完整口径见 `reports/04_asymmetric_mapping_results.md`。

### 8. 同一结构–电子化学空间

对同时具有正带隙、电子/空穴有效质量和介电响应的 3,259 个 JARVIS 半导体，现已构造
“电子感知结构空间”：结构 superblock 为 162 维组成 + 147 维纯几何 SOAP，电子 superblock
为四个电子性质的 OOF 预测；两个 superblock 分别归一到相同总惯量后 1:1 融合。

联合二维空间保留 51.3% 方差。joint axis 1 主要体现电子梯度（`Eg` +0.70、`log ε` −0.77、
平均 `log m*` +0.55），joint axis 2 主要体现几何梯度（SOAP PC1 +0.95、`log V/atom` +0.75）。
这说明同一空间可以同时表达两类信息，但点云是连续流形，不是两个天然簇。图和读图说明见
`reports/05_unified_chemical_space.md`。图中另以青色星形在左右两个镜头同时标出 14 个经典热电
基准化学式（27 个结构条目）；它们只是外部参考地标，不是模型给出的候选排名。

### 9. 结构邻域与电子结构邻域的交集

为避免“同一张图只有梯度、没有交集”，现以 14 个经典热电化学式为参考，分别在
`composition + SOAP` 结构/化学空间和 `Eg + m*e + m*h + ε` 电子结构空间计算到最近三个
参考家族质心的平均距离。排除参考条目后，各取最相似的前 5%。

在 3,232 个非参考材料中，两边各有 162 个点，实际交集为 19 条结构（18 个化学式），高于
随机期望 8.12，富集 2.34×（超几何 `p=3.57×10⁻⁴`）；逐一删除一个参考家族后，14 条仍以
至少 80% 的频率留在交集中。交集以 Bi/Sb/Tl–Te/Se 硫属化物为主。它表示“双描述符空间都
接近已知热电基准”，不是 `zT` 或独立 `κL` 排名。方法、图和完整口径见
`reports/06_te_reference_intersection.md`。

## 关键审计修正

1. **共同队列 kNN**：每个视图对先取共同材料，再在共同队列内重建邻域；旧实现先在全集找邻居
   再丢掉非共同邻居，使 Eg 的随机基线错误。
2. **稳定并列破平局**：ID 次序键改为 BLAKE2b 64-bit；旧 little-endian 前四字节对 `mp-*`
   大面积碰撞。打乱行序后的最大偏差现在低于 `2e-16`。
3. **实验目标不复制**：4,401 个实验规范化化学式中，仅 59 个能唯一映射到单一 MP material_id；
   另有 44 个歧义化学式。材料级消融只用这 59 条。
4. **标定不可跨生成器倒推**：Eg 只能使用 Eg-1D 专属模拟曲线；其逆映射范围约
   0.086–0.173 只是该模拟器参数，不是“物理等效 R²”。旧版 0.48 已撤回。
5. **增量检验修正**：块严格分离、RF 使用全部特征、偏 Spearman 在秩上带截距残差化；
   原“条件互信息”改名为交叉拟合残差 MI 代理并报告置换偏置，不能解释为严格 CMI。
6. **筛选 fail-closed**：五个必要条件不满足时不再生成候选排名。

## 数据口径

| 数据 | 口径 | 当前用途 |
|---|---|---|
| MP elasticity | 12,246→12,156（90 个唯一异常材料剔除） | 结构、弹性、Clarke/Cahill 下界、Snyder 模型 |
| MP summary | 2,472 个共同材料有可用 Eg 视图 | material_id 级电子带隙对照 |
| JARVIS dft_3d | 7,104 个唯一化学式有 600 K CRTA 输运 | κe/PF/σ/S 的探索性公式层分析 |
| Starrydata2 | 4,401 个聚合化学式；96 个描述符交集；70 个输运交集 | 实验 κL 验证 |

JARVIS 的 σ/τ 与 κe/τ 数值尺度未转成绝对 zT 分母；这里只使用秩和对数距离。方法背景见
[JARVIS-DFT 数据论文](https://www.nature.com/articles/s41524-020-00440-1)、
[JARVIS 数据字段说明](https://jarvis-materials-design.github.io/dbdocs/jarvisdft/)和
[pymatgen BoltzTraP 文档](https://pymatgen.org/pymatgen.electronic_structure.html)。

## 主要输出

- `mp_kappaL/processed/channel_property_correlations.csv`：各性质对 κe/模型 κL/实验 κL 的相关与 CI
- `mp_kappaL/processed/channel_space_distance.csv`：共同性质空间的距离相关
- `mp_kappaL/processed/channel_space_overlap.csv`：共同队列 10-NN 重叠、随机基线与富集
- `mp_kappaL/processed/ke_kl_crosschannel.csv`：κe–κL 直接交叉验证
- `mp_kappaL/processed/experimental_subset_sensitivity.csv`：多晶型/实验曲线子集与留一法
- `mp_kappaL/processed/geometry_increment_summary.csv`：严格配对 ΔR²
- `mp_kappaL/processed/screening_readiness.csv`：筛选阻断条件
- `mp_kappaL/figures/channel_property_correlations.png`
- `mp_kappaL/figures/channel_space_geometry.png`
- `mp_kappaL/figures/channel_space_pca.png`
- `mp_kappaL/figures/ke_kl_crosschannel.png`
- `mp_kappaL/figures/corrected_audit_summary.png`
- `mp_kappaL/figures/structure_space_channel_clusters.png`：**历史诊断图，不作材料结论**
- `mp_kappaL/figures/structure_space_two_clusters.png`：**历史诊断图，标签口径不适合筛选**
- `mp_kappaL/processed/structure_channel_membership.csv`
- `mp_kappaL/processed/structure_channel_separation.csv`
- `mp_kappaL/processed/structure_channel_pair_overlap.csv`
- `reports/03_asymmetric_structure_to_transport_plan.md`：不对称映射、品质因子和泄漏审计的修订方案
- `mp_kappaL/figures/asymmetric_mapping_decisive_test.png`：分组 CV、几何增量与泄漏对照主图
- `mp_kappaL/processed/asymmetric_mapping_summary.csv`：OOF 汇总分数与 group-bootstrap CI
- `mp_kappaL/processed/asymmetric_mapping_folds.csv`：逐折分数
- `mp_kappaL/processed/asymmetric_mapping_oof.parquet`：逐材料 OOF 预测/残差
- `reports/04_asymmetric_mapping_results.md`：判决性实验结论
- `mp_kappaL/figures/unified_structure_electronic_space.png`：同一坐标的结构/电子双镜头图
- `mp_kappaL/processed/unified_chemical_space_coordinates.csv`：逐材料联合坐标与性质
- `mp_kappaL/processed/unified_chemical_space_te_benchmarks.csv`：图中基准热电体系及代表点标记
- `mp_kappaL/processed/unified_chemical_space_weight_sensitivity.csv`：块权重敏感性
- `reports/05_unified_chemical_space.md`：构造方法、读图方式与限制
- `mp_kappaL/figures/te_reference_dual_space_intersection.png`：结构/电子参考邻域的并集与交集
- `mp_kappaL/processed/te_reference_dual_space_membership.csv`：逐材料双空间相似度与交集稳定性
- `reports/06_te_reference_intersection.md`：交集定义、统计结果与正确解释

## 复现

环境为 Python 3.11 的 `te_manifold`；运行前把仓库根目录加入 `PYTHONPATH`。

```bash
export PYTHONPATH=$PWD

# 21 个回归/循环性/对齐测试
python -m pytest tests -q

# 基础视图和目标
python -m mp_kappaL.build_views
python -m mp_kappaL.build_electronic
python -m mp_kappaL.kappa_L_targets

# 跨视图、标定和增量信息
python -m mp_kappaL.crossview_analysis
python -m mp_kappaL.row_order_audit
python -m mp_kappaL.calibration
python -m mp_kappaL.ablation
python -m mp_kappaL.soap_robustness

# κe–κL 空间关系、筛选阻断和总览图
python -m mp_kappaL.dual_channel
python -m mp_kappaL.channel_space_analysis
# 历史诊断图；默认复现不再运行，不能用于候选结论：
# python -m mp_kappaL.structure_channel_map
python -m mp_kappaL.screening
python -m mp_kappaL.audit_summary

# 不对称结构→电子性质判决性实验（默认读取已生成的 JARVIS SOAP 缓存）
python -m mp_kappaL.asymmetric_mapping

# 等权结构–电子联合化学空间
python -m mp_kappaL.unified_chemical_space

# 已知热电参考家族的结构/电子双空间交集
python -m mp_kappaL.te_reference_intersection
```

`crossview_analysis` 是最耗时步骤；需要数 GB 内存。`channel_space_analysis` 固定抽取 2,500 个
化学式做距离矩阵，避免 O(N²) 内存失控。所有抽样、bootstrap 和并列处理均固定种子。

更详细的审计记录见 `reports/audit_corrected.md`，候选阻断状态见 `reports/candidates.md`。
