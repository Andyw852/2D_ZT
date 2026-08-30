# JARVIS 二维热电输运图谱 Phase L0：二维输运主值谱物理一致性审计与最终特征冻结

## 一、任务定位

你现在是我的二维热电材料、材料信息学、Boltzmann 输运、张量分析、统计学习和流形学习研究助手。

前两阶段的数据获取与输运数据审计已经完成。

当前不要立即开始 SOAP、UMAP、Diffusion Map、图嵌入或多流形融合。

在进入正式结构流形建模之前，需要额外完成一个关键阶段：

**Phase L0：二维输运主值谱的物理一致性审计与最终 Transport Representation 冻结。**

本阶段的目标不是得到漂亮的降维图，而是解决一个更基础的问题：

> JARVIS 给出的 Seebeck、电导率、电子热导率和有效质量均为 3×3 张量对角化后的 3 个本征值，但没有本征向量，因此我们究竟应该怎样把这 3 个主值转化为适合二维材料比较的、排列不变且具有明确物理含义的特征？

本阶段完成后，需要明确冻结：

- n-type transport feature set；
- p-type transport feature set；
- electronic feature set；
- PF 的使用方式；
- κ_e 的使用方式；
- 各向异性描述符；
- Seebeck 符号异常的处理原则。

完成 Phase L0 后立即 STOP。

不要进入真正的流形建模。

---

# 二、当前已经确认的数据事实

当前数据库：

`JARVIS dft_2d`

当前材料总数：

`1103`

通过 JARVIS 官方 NIST 静态 XML：

`https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/{JID}`

成功获得全部 1103 个材料记录。

没有下载失败。

---

## 2.1 已恢复的输运数据

成功恢复：

`807` 条具有原始输运主值谱的数据。

其中：

- n-type 完整 S / σ / κ_e 集合：806 条；
- p-type 完整 S / σ / κ_e 集合：803 条；
- effective mass：678 条。

---

## 2.2 输运条件

已经通过 jarvis-tools 源码确认：

所有这些输运数据都来自：

`T = 600 K`

以及：

`|carrier concentration| = 1e20 cm^-3`

其中：

- n-type 对应电子掺杂；
- p-type 对应空穴掺杂。

这些数据不是：

- 多温度平均；
- doping 扫描平均；
- transport response surface。

而是固定条件下的输运张量结果。

---

## 2.3 原始张量信息层级

JARVIS XML 中保存的是：

`3 eigenvalues`

而不是完整：

`3 × 3 tensor`

也没有：

`eigenvectors`

因此当前知道的是：

**principal-value spectrum**

但不知道：

**principal directions**

因此当前不能把 3 个主值解释为：

`x / y / z`

也不能解释为：

`xx / yy / zz`

---

# 三、本阶段最重要的物理原则

## 原则 1：三个本征值的数组顺序没有可靠物理意义

禁止直接把：

`[eig1, eig2, eig3]`

理解成：

`[x, y, z]`

或者：

`[in-plane-1, in-plane-2, out-of-plane]`

因为 XML 中不存在 eigenvectors。

---

## 原则 2：特征必须尽量 permutation invariant

所谓 permutation invariant，是指：

交换三个本征值顺序以后，最终材料描述符不应该发生变化。

例如可以使用：

- mean；
- median；
- std；
- MAD；
- min；
- max；
- range；
- principal-value ratios。

而不是直接依赖：

`eig1`

`eig2`

`eig3`

的原始排列。

---

## 原则 3：二维材料不能简单把最小本征值称为 z 方向

二维材料通常存在：

- 两个较强的面内输运通道；
- 一个较弱的面外输运通道。

但仅仅因为：

`λ1 < λ2 < λ3`

不能直接得出：

`λ1 = out-of-plane`

因为没有 eigenvectors。

因此后续只能使用更加谨慎的名称：

- weakest principal channel；
- middle principal channel；
- strongest principal channel；
- dominant-channel pair；
- suppressed-channel contrast。

禁止直接写：

`out-of-plane conductivity`

除非以后获得原始 tensor 或 eigenvectors。

---

## 原则 4：PF 默认是外部性能标签，而不是主流形输入

PF 来源于：

`PF = S² × σ / 1e6`

因此它不是独立基础输运自由度。

主 Transport View 默认由：

`Seebeck + conductivity`

以及必要的谱各向异性 descriptor 构成。

PF 主要作为：

**external performance label**

用于验证：

> 基础输运空间中是否自然出现 high-PF 区域。

---

## 原则 5：κ_e 默认不和 σ 同时等权进入主模型

当前已经发现：

n-type：

`Pearson(log σ, log κ_e) ≈ 0.960`

`Spearman(log σ, log κ_e) ≈ 0.928`

p-type：

`Pearson(log σ, log κ_e) ≈ 0.965`

`Spearman(log σ, log κ_e) ≈ 0.947`

说明二者高度相关。

因此：

主模型优先使用：

`S + σ`

κ_e 作为 sensitivity model。

禁止因为名字不同，就把 σ 和 κ_e 当作两个完全独立物理层等权计入。

---

# 四、Phase L0 总体流程

严格按照以下顺序执行：

```text
L0-A
PF 主值配对歧义审计
        ↓
L0-B
σ 和 κ_e 的主值谱二维性审计
        ↓
L0-C
Seebeck 主值谱与符号一致性审计
        ↓
L0-D
mean / full-spectrum / dominant-spectrum 邻域比较
        ↓
L0-E
effective-mass 主值谱审计
        ↓
L0-F
Transport representation 冗余比较
        ↓
L0-G
冻结最终 n/p Transport Features
        ↓
L0-H
冻结 Electronic Features
        ↓
STOP
```

每一个子阶段完成后都要：

1. 执行程序；
2. 检查异常；
3. 保存 CSV / Parquet；
4. 保存图；
5. 写 Markdown 小结；
6. 再进入下一步。

---

# 五、L0-A：PF Principal-Value Pairing Ambiguity

这是本阶段优先级最高的问题之一。

JARVIS 源码逻辑是：

```python
S_eigs = eigvals(S_tensor)
sigma_eigs = eigvals(sigma_tensor)

PF_eigs = S_eigs**2 * sigma_eigs / 1e6
```

但需要注意：

`np.linalg.eigvals()`

分别作用于 S tensor 和 σ tensor。

除非两个 tensor：

- 具有共同本征向量；
- 且 eigvals 返回顺序保持对应；

否则：

`S_eigs[i]`

和：

`sigma_eigs[i]`

未必属于同一个物理 principal direction。

因此需要定量评估：

**PF 对 eigenvalue pairing 顺序有多敏感。**

---

# 六、L0-A1：生成 6 种所有可能配对

对于每个材料：

Seebeck 三主值：

`S = [S1, S2, S3]`

conductivity 三主值：

`sigma = [σ1, σ2, σ3]`

固定 S 顺序。

对 σ 的所有：

`3! = 6`

种 permutation 进行枚举。

对于 permutation p：

分别计算：

`PF_i(p) = S_i² × σ_p(i) / 1e6`

然后：

`PF_mean(p) = mean(PF_1, PF_2, PF_3)`

获得：

`PF_mean_perm_1`

到：

`PF_mean_perm_6`

---

# 七、L0-A2：定义 PF 配对歧义指标

对于每个材料计算：

`PF_pair_min = min(PF_mean_perm_1 ... PF_mean_perm_6)`

`PF_pair_max = max(...)`

`PF_pair_median = median(...)`

定义：

`PF_pairing_ambiguity = (PF_pair_max - PF_pair_min) / (|PF_pair_median| + epsilon)`

其中 epsilon 不要随意写死。

优先使用：

`epsilon = 全部非零 PF_pair_median 的最小绝对值 / 10`

并进行敏感性检查。

---

# 八、L0-A3：统计 PF pairing ambiguity

分别对：

- n-type；
- p-type；

统计：

- median；
- mean；
- 90th percentile；
- 95th percentile；
- max。

同时统计：

`ambiguity < 0.05`

`0.05 ≤ ambiguity < 0.20`

`0.20 ≤ ambiguity < 0.50`

`ambiguity ≥ 0.50`

的材料比例。

---

# 九、L0-A4：与 JARVIS PF 比较

计算 JARVIS 原始 PF spectrum 和：

6 种 permutation 中每一种结果的差异。

判断 JARVIS 当前 pairing 是否：

- 位于所有可能 pairing 的中间；
- 经常接近极大值；
- 经常接近极小值；
- 对 majority materials 几乎无影响。

---

# 十、L0-A5：PF 使用决策

根据结果进行分类。

如果：

大多数材料：

`PF_pairing_ambiguity < 0.05`

则：

PF pairing ambiguity 可以认为较小。

PF 可以较放心作为 external performance label。

如果：

大量材料：

`PF_pairing_ambiguity > 0.20`

则：

必须明确将 JARVIS PF 定义为：

**database-defined PF under JARVIS eigenvalue pairing convention**

而不能解释成严格 crystallographic principal-direction PF。

PF 仍可用于：

- 数据库内部排序；
- high-PF mapping；
- relative comparison。

但不能过度解释方向性。

---

# 十一、L0-A6：输出

生成：

`data/audit/pf_pairing_ambiguity_n.csv`

`data/audit/pf_pairing_ambiguity_p.csv`

`reports/pf_pairing_ambiguity.md`

图：

`figures/pf_pairing_ambiguity_hist_n.png`

`figures/pf_pairing_ambiguity_hist_p.png`

`figures/pf_pairing_vs_pf_n.png`

`figures/pf_pairing_vs_pf_p.png`

---

# 十二、L0-B：Conductivity Principal Spectrum 审计

对于每个材料的 conductivity eigenvalues：

先确保：

`σ_i > 0`

或者记录异常。

然后排序：

`σ_(1) ≤ σ_(2) ≤ σ_(3)`

注意：

这里的 1、2、3 只表示数值排序。

不是 x/y/z 方向。

---

# 十三、L0-B1：定义 Conductivity Scale

不要直接用三个 raw values。

定义：

### Arithmetic mean

`sigma_mean`

### Median

`sigma_median`

### Geometric mean of all channels

`sigma_geo_all = (σ1 × σ2 × σ3)^(1/3)`

只有全部严格正时计算。

### Dominant-channel geometric scale

`sigma_dom_geo = sqrt(σ2 × σ3)`

这个量代表：

两个较强 principal channels 的整体输运尺度。

但命名必须保持：

**dominant-channel conductivity scale**

禁止直接叫：

`in-plane conductivity`

---

# 十四、L0-B2：定义谱各向异性

定义：

### Overall spectral anisotropy

`A_sigma_total = log10(σ3 / σ1)`

### Suppressed-channel contrast

`D_sigma = log10(σ2 / σ1)`

### Dominant-channel anisotropy

`A_sigma_dom = log10(σ3 / σ2)`

这样：

`A_sigma_total ≈ D_sigma + A_sigma_dom`

其中：

`D_sigma`

更可能包含 quasi-2D dimensionality 信息。

而：

`A_sigma_dom`

更接近两个主要传输通道之间的差异。

注意：

不能直接称：

`D_sigma = out-of-plane anisotropy`

也不能严格称：

`A_sigma_dom = in-plane anisotropy`

更安全的术语：

- suppressed-channel contrast；
- dominant-principal-channel anisotropy。

---

# 十五、L0-B3：统计二维性模式

分析全体二维材料中：

`σ1 / σ2`

以及：

`σ2 / σ3`

的分布。

重点回答：

1. 是否大量材料存在：

`σ1 ≪ σ2 ≈ σ3`

？

2. 是否存在：

`σ1 ≈ σ2 ≈ σ3`

的近似 isotropic materials？

3. 是否存在：

`σ1 ≪ σ2 ≪ σ3`

的高度各向异性 materials？

对材料进行谱模式分类：

### Type C1

`σ1 ≪ σ2 ≈ σ3`

### Type C2

`σ1 < σ2 < σ3`

### Type C3

`σ1 ≈ σ2 ≈ σ3`

### Type C4

其他异常模式。

阈值不能直接拍脑袋。

先基于：

`log10(σ2/σ1)`

和：

`log10(σ3/σ2)`

的总体分布决定。

---

# 十六、L0-B4：κ_e 做相同处理

对于：

`κe_(1) ≤ κe_(2) ≤ κe_(3)`

定义：

`kappa_mean`

`kappa_median`

`kappa_dom_geo`

`D_kappa`

`A_kappa_dom`

并比较：

`D_sigma vs D_kappa`

`A_sigma_dom vs A_kappa_dom`

计算：

- Pearson；
- Spearman。

如果二者几乎完全相关：

进一步证明 κ_e 可以作为 sensitivity variable 而不是独立主自由度。

---

# 十七、L0-B5：输出

保存：

`data/audit/conductivity_spectrum_audit_n.csv`

`data/audit/conductivity_spectrum_audit_p.csv`

`data/audit/kappa_spectrum_audit_n.csv`

`data/audit/kappa_spectrum_audit_p.csv`

图：

`figures/sigma_dimensionality_n.png`

`figures/sigma_dimensionality_p.png`

`figures/sigma_D_vs_A_n.png`

`figures/sigma_D_vs_A_p.png`

`figures/kappa_D_vs_A_n.png`

`figures/kappa_D_vs_A_p.png`

---

# 十八、L0-C：Seebeck Principal Spectrum 审计

Seebeck 与 conductivity 不同。

Seebeck 可以：

- 为正；
- 为负；
- 接近零；
- 三个 principal values 符号不一致。

因此不能使用：

`max / min`

作为 anisotropy。

---

# 十九、L0-C1：Seebeck 基础谱特征

对三个本征值：

`S1, S2, S3`

构造排列不变的：

`S_mean`

`S_median`

`S_std`

`S_MAD`

`S_min`

`S_max`

`S_range`

`S_abs_mean`

`S_abs_max`

定义：

`S_relative_spread = S_std / (|S_mean| + epsilon_S)`

同时考虑：

`S_MAD`

作为对单个异常 eigenvalue 更稳健的 dispersion measure。

---

# 二十、L0-C2：定义 Sign Consistency

对 n-type：

理论常见符号：

`S < 0`

计算：

`N_expected_sign_n`

其值可能为：

`0, 1, 2, 3`

分别表示三个 principal values 中有多少个为负。

定义：

`sign_fraction_n = N_expected_sign_n / 3`

p-type 同样：

期望：

`S > 0`

定义：

`N_expected_sign_p`

以及：

`sign_fraction_p`

---

# 二十一、L0-C3：按符号一致性分类

n-type：

### N3

3/3 为负。

### N2

2/3 为负。

### N1

1/3 为负。

### N0

0/3 为负。

p-type：

### P3

3/3 为正。

### P2

2/3 为正。

### P1

1/3 为正。

### P0

0/3 为正。

分别统计材料数量和比例。

---

# 二十二、L0-C4：不要直接把符号异常归因于面外分量

当前没有 eigenvectors。

因此禁止直接写：

> 第三个符号异常主值是 out-of-plane contribution。

这是没有证据的。

需要验证至少三个可能解释：

### Hypothesis H1

suppressed-channel numerical instability。

### Hypothesis H2

small-gap / bipolar transport。

### Hypothesis H3

multiband complexity。

当前数据不一定足以完全证明三者。

但至少可以寻找统计证据。

---

# 二十三、L0-C5：Sign Consistency 与 Band Gap

对：

`N3/N2/N1/N0`

以及：

`P3/P2/P1/P0`

分别统计：

`OptB88vdW band gap`

的：

- mean；
- median；
- IQR；
- distribution。

画：

`sign consistency vs band gap`

重点检查：

符号严重违例：

`N0/N1`

或：

`P0/P1`

是否显著集中于：

`small band gap`

区域。

---

# 二十四、L0-C6：Sign Consistency 与其他量

进一步比较：

- PF；
- effective mass；
- sigma dimensionality contrast；
- sigma dominant anisotropy；
- chemical family。

目标：

判断符号异常更可能与：

- 低带隙；
- 多带输运；
- 极弱 transport channel；

哪一种因素相关。

注意：

只能说：

`associated with`

不能在没有充分证据时写：

`caused by`

---

# 二十五、L0-C7：比较 mean 和 median

这是非常重要的。

比较：

`S_mean`

和：

`S_median`

对于 sign-violation materials 的稳定性。

计算：

`sign(S_mean)`

`sign(S_median)`

与 majority eigenvalue sign 的一致程度。

定义：

`majority_sign = sign of at least 2 of 3 eigenvalues`

统计：

### mean accuracy

`sign(S_mean) == majority_sign`

的比例。

### median accuracy

`sign(S_median) == majority_sign`

的比例。

因为 median 对一个异常 eigenvalue 天然稳健。

---

# 二十六、L0-C8：Seebeck representation 候选

建立三种候选。

### S-Model 1：Mean only

`[S_mean]`

### S-Model 2：Median based

`[S_median]`

### S-Model 3：Spectrum robust

`[S_median, S_MAD, S_range, sign_fraction]`

必要时保留：

`S_mean`

作为辅助变量。

后续比较它们对材料邻域的影响。

---

# 二十七、L0-C9：输出

生成：

`data/audit/seebeck_spectrum_n.csv`

`data/audit/seebeck_spectrum_p.csv`

`data/audit/seebeck_sign_consistency_n.csv`

`data/audit/seebeck_sign_consistency_p.csv`

`reports/seebeck_sign_analysis.md`

图：

`figures/seebeck_sign_vs_gap_n.png`

`figures/seebeck_sign_vs_gap_p.png`

`figures/seebeck_mean_vs_median_n.png`

`figures/seebeck_mean_vs_median_p.png`

---

# 二十八、L0-D：比较三种 Transport Spectrum Representation

现在需要回答：

> 之前发现 full spectrum 会显著改变材料邻域，这究竟来自真正有价值的 dominant-channel anisotropy，还是主要因为 weakest channel 与两个主要通道差别巨大？

因此建立三个模型。

---

# 二十九、Representation R1：Mean-Only

n-type：

`[S_mean, log_sigma_mean]`

p-type：

`[S_mean, log_sigma_mean]`

这是 OPTIMADE scalar representation 的近似。

---

# 三十、Representation R2：Full Spectrum Statistical

使用：

`S_median`

`S_MAD`

`S_range`

`log_sigma_mean`

`D_sigma`

`A_sigma_dom`

这是不依赖原始 eigenvalue 顺序的完整 spectrum descriptor。

---

# 三十一、Representation R3：Dominant-Channel Focused

使用：

`S_median`

`S_MAD`

`log_sigma_dom_geo`

`A_sigma_dom`

同时保留：

`D_sigma`

但作为 dimensionality descriptor。

这一版本的目标是：

减少 weakest channel 对平均量的强烈污染。

注意：

`dominant channels`

不能直接称：

`in-plane channels`

因为没有 eigenvectors。

---

# 三十二、L0-D1：标准化规则

Seebeck-like features：

使用 StandardScaler 或 RobustScaler。

sigma scale：

先：

`log10(value)`

再标准化。

anisotropy log ratios：

本身已经是 logarithmic descriptor，

不要再次 log。

所有 scaler 独立保存。

---

# 三十三、L0-D2：比较材料距离

分别基于：

R1、R2、R3

计算：

pairwise distance。

不要做最终 UMAP。

这里只进行 representation audit。

计算：

- Spearman rank correlation of distances；
- Pearson correlation of distances。

---

# 三十四、L0-D3：比较 kNN

测试：

`k = 5, 10, 20, 30`

分别计算：

`R1 vs R2`

`R1 vs R3`

`R2 vs R3`

的：

kNN overlap。

重点看：

R2 和 R3 是否高度一致。

如果：

`overlap(R2,R3) > 0.8`

而：

`overlap(R1,R2) ≈ 0.2`

说明：

主要问题确实来自 mean-only 丢失 spectrum 信息。

如果：

R2 和 R3 差异仍然很大，

说明 weakest channel 本身在材料分类中具有非常强的影响。

---

# 三十五、L0-D4：PF Smoothness Preview

虽然这一阶段不做正式 graph manifold，

但允许用简单 kNN neighborhood 做 preliminary PF smoothness comparison。

分别使用：

R1、R2、R3

检查：

high PF materials 是否拥有更相似的局部邻域。

不要把结果作为最终论文结论。

这里只用于选择 transport representation。

---

# 三十六、L0-D5：最终 representation 选择原则

优先选择：

1. permutation invariant；
2. 对单个异常主值稳健；
3. 保留显著 anisotropy；
4. 不被 weakest channel 完全支配；
5. 与 PF landscape 有合理关联；
6. n/p 两套处理逻辑一致；
7. feature 数量尽量简洁。

不要因为 R2 特征最多就自动选择 R2。

---

# 三十七、L0-E：Effective-Mass Spectrum Audit

当前 effective mass：

`N = 678`

同样是三个 principal values。

首先检查字段究竟对应：

- electron effective mass；
- hole effective mass；
- 或其他数据库定义。

以实际 XML 字段为准。

如果 n/p 有独立记录：

必须分别处理。

---

# 三十八、L0-E1：Effective-Mass 数值审计

检查：

- positive；
- negative；
- zero；
- Inf；
- NaN；
- extreme values。

如果存在负 effective mass：

先检查 JARVIS 定义。

不要自动 abs。

必须确认：

它是否已经是 carrier effective mass magnitude，

还是 band-curvature signed quantity。

---

# 三十九、L0-E2：Permutation-Invariant Effective-Mass Features

如果物理定义允许使用绝对值：

可构造：

`m_abs_mean`

`m_abs_median`

`m_abs_std`

`m_abs_min`

`m_abs_max`

`m_abs_range`

以及：

`m_spectral_ratio = log10(m_abs_max / m_abs_min)`

如果不能使用绝对值：

需要根据实际物理定义另外设计。

禁止未经确认直接处理。

---

# 四十、L0-E3：Effective Mass 与 Transport 对应性

如果存在 electron effective mass：

与：

n-Transport

比较。

如果存在 hole effective mass：

与：

p-Transport

比较。

分别计算：

- correlation with S；
- correlation with sigma；
- correlation with PF；
- neighborhood consistency preview。

目的：

确定 effective mass 是否应该进入：

Electronic View

而不是 Transport View。

默认优先：

**effective mass 属于 Electronic View。**

---

# 四十一、L0-F：Transport Feature 冗余与敏感性模型

经过前面的 spectrum audit 后，

建立三个 Transport candidates。

---

# 四十二、Transport Model T1：主模型

默认：

`robust Seebeck spectrum`

+

`conductivity spectrum`

例如候选：

`S_median`

`S_MAD`

`S_range`

`log_sigma_dom_geo`

`D_sigma`

`A_sigma_dom`

不输入：

`PF`

不输入：

`κ_e`

---

# 四十三、Transport Model T2：κ_e 替代模型

把 conductivity 部分替换为：

`κ_e spectrum`

例如：

`S_median`

`S_MAD`

`S_range`

`log_kappa_dom_geo`

`D_kappa`

`A_kappa_dom`

用来判断：

使用 σ 或 κ_e 是否实际上得到几乎相同的材料空间。

---

# 四十四、Transport Model T3：σ + κ_e 全部加入

只用于 sensitivity analysis：

`Seebeck descriptors`

+

`conductivity descriptors`

+

`κ_e descriptors`

不要作为默认模型。

---

# 四十五、T1 / T2 / T3 比较

计算：

- pairwise-distance rank correlation；
- kNN overlap；
- PCA effective rank；
- PF neighborhood smoothness preview。

如果：

`T1 vs T3 kNN overlap > 0.9`

说明加入 κ_e 几乎不改变材料邻域。

则正式主模型使用：

`T1`

更加简洁。

如果：

T1 和 T2 也高度相似，

则进一步确认：

σ 和 κ_e 在当前固定条件下属于同一个主要 transport degree of freedom。

---

# 四十六、PF 永远不加入 T1 主模型

PF 继续作为：

`external performance label`

后续真正建立流形后，

分别映射：

`PF_n`

和：

`PF_p`

检验：

> high-PF materials 是否自然集中于某些 transport / joint manifold 区域。

---

# 四十七、L0-G：冻结最终 n-Type Transport Features

完成所有审计以后，

生成：

`features/transport/n_transport_features_v1.parquet`

文件只保留最终正式使用的 n-type transport features。

同时生成：

`features/transport/n_transport_features_candidates.parquet`

保留 T1/T2/T3 的候选变量。

---

# 四十八、L0-G1：n-Type Metadata

为每个 feature 建立：

`feature_name`

`physical_meaning`

`source_property`

`transform`

`unit`

`permutation_invariant`

`used_in_main_model`

`reason`

保存：

`features/transport/n_transport_feature_metadata.csv`

---

# 四十九、L0-G2：冻结 p-Type Transport Features

完全按照相同原则生成：

`features/transport/p_transport_features_v1.parquet`

以及：

`features/transport/p_transport_feature_metadata.csv`

n/p 两套结构尽量保持对称。

如果由于数据异常导致变量不同，

必须在报告中解释。

---

# 五十、L0-H：冻结 Electronic View

Electronic View 不要塞入所有输运量。

其目标是描述：

**band/electronic structure state**

而不是直接描述 transport performance。

第一版候选：

`OptB88vdW band gap`

+

`effective-mass spectral descriptors`

如果 effective mass 具有 n/p 分离信息：

可分别建立：

`Electronic-n`

和：

`Electronic-p`

候选 representation。

---

# 五十一、Band Gap 的角色

OptB88vdW band gap：

覆盖率：

`100%`

因此可以作为 Electronic View 的基础变量。

MBJ gap：

覆盖率：

约 22%。

暂时不要进入主 Electronic View。

将其作为：

**higher-level electronic validation layer**

HSE06 gap：

覆盖率太低。

仅 exploratory。

---

# 五十二、不要让 Electronic View 被 Eg 一个变量主导

如果 effective mass 只覆盖 678 个材料，

不要简单建立：

`Electronic = [Eg, m*]`

然后把其余 425 个材料全部删除。

需要保存两层概念：

### Eg electronic layer

覆盖 1103。

### Rich electronic layer

`Eg + effective-mass spectrum`

覆盖 678。

下一阶段做 partial multi-view graph 时决定如何组合。

当前只准备 features。

---

# 五十三、最终推荐四大物理 View

本阶段结束时，目标架构优先为：

```text
                         Structure
                            |
          +-----------------+-----------------+
          |                 |                 |
      Electronic       n-Transport       p-Transport
```

其中：

### Structure

未来：

SOAP + composition + local geometry + symmetry。

### Electronic

Eg + effective-mass spectrum。

### n-Transport

robust Seebeck spectrum + conductivity spectrum。

### p-Transport

robust Seebeck spectrum + conductivity spectrum。

PF：

external performance label。

κ_e：

sensitivity representation。

---

# 五十四、本阶段不要把这四个 View 融合

本阶段只负责：

**定义和冻结 representation。**

禁止开始：

- kNN graph construction；
- UMAP；
- Diffusion Map；
- Laplacian Eigenmap；
- manifold alignment；
- supra adjacency；
- joint embedding。

这些属于下一阶段。

---

# 五十五、本阶段必须生成的核心输出

至少生成：

## PF pairing

`data/audit/pf_pairing_ambiguity_n.csv`

`data/audit/pf_pairing_ambiguity_p.csv`

## Conductivity / kappa spectrum

`data/audit/conductivity_spectrum_audit_n.csv`

`data/audit/conductivity_spectrum_audit_p.csv`

`data/audit/kappa_spectrum_audit_n.csv`

`data/audit/kappa_spectrum_audit_p.csv`

## Seebeck

`data/audit/seebeck_spectrum_n.csv`

`data/audit/seebeck_spectrum_p.csv`

`data/audit/seebeck_sign_consistency_n.csv`

`data/audit/seebeck_sign_consistency_p.csv`

## Effective mass

`data/audit/effective_mass_spectrum.csv`

## Representation comparison

`data/audit/transport_representation_comparison_n.csv`

`data/audit/transport_representation_comparison_p.csv`

`data/audit/transport_T1_T2_T3_comparison.csv`

---

# 五十六、最终 Feature 文件

必须生成：

`features/transport/n_transport_features_v1.parquet`

`features/transport/p_transport_features_v1.parquet`

`features/transport/n_transport_feature_metadata.csv`

`features/transport/p_transport_feature_metadata.csv`

以及：

`features/electronic/electronic_features_v1.parquet`

`features/electronic/electronic_feature_metadata.csv`

---

# 五十七、报告文件

至少输出：

`reports/pf_pairing_ambiguity.md`

`reports/seebeck_sign_analysis.md`

`reports/transport_spectrum_physics.md`

`reports/transport_representation_selection.md`

`reports/electronic_feature_selection.md`

`reports/phase_L0_summary.md`

---

# 五十八、建议绘制的图

至少生成：

1. PF pairing ambiguity distribution，n-type。

2. PF pairing ambiguity distribution，p-type。

3. `D_sigma vs A_sigma_dom`，n-type。

4. `D_sigma vs A_sigma_dom`，p-type。

5. Seebeck sign consistency vs band gap，n-type。

6. Seebeck sign consistency vs band gap，p-type。

7. `S_mean vs S_median`。

8. R1 / R2 / R3 kNN overlap。

9. T1 / T2 / T3 neighbor overlap。

10. Effective-mass spectral descriptor distribution。

图必须：

- 白底；
- 标签清楚；
- 单位完整；
- 不使用未经验证的 x/y/z 方向名称。

---

# 五十九、Phase L0 最终必须回答的问题

完成以后，必须明确回答以下问题。

1. JARVIS PF 对 S 与 σ 的 eigenvalue pairing 是否敏感？

2. PF 是否可以继续作为可靠的数据库内部 external performance label？

3. conductivity spectrum 是否普遍表现为：

`weakest channel << two dominant channels`

？

4. 当前 full-spectrum anisotropy 主要反映：

suppressed-channel contrast，

还是 dominant-channel anisotropy？

5. 两个 dominant conductivity channels 之间的差异有多大？

6. κ_e 与 σ 的 spectrum anisotropy 是否也高度一致？

7. 15% 左右的 Seebeck sign violations 是否与 small band gap 显著相关？

8. `S_median` 是否比 `S_mean` 更符合 majority principal-value sign？

9. mean-only、full-spectrum、dominant-spectrum 三种 representation 的材料邻域差别是多少？

10. anisotropy spectrum 信息是否应该保留？

11. T1、T2、T3 哪个最适合作为主 Transport representation？

12. κ_e 是否只需要作为 sensitivity variable？

13. PF 是否继续不进入主 Transport View？

14. effective-mass spectrum 最终应如何表示？

15. 最终 Electronic View 使用哪些变量？

16. 最终 n-Transport View 使用哪些变量？

17. 最终 p-Transport View 使用哪些变量？

---

# 六十、必须给出最终 Feature Decision Table

生成类似：

| View | Feature | Keep in main model? | Role | Reason |
|---|---|---|---|---|
| n-Transport | S_median | ? | primary | ? |
| n-Transport | S_MAD | ? | anisotropy | ? |
| n-Transport | S_range | ? | spectrum spread | ? |
| n-Transport | log_sigma_dom_geo | ? | transport scale | ? |
| n-Transport | D_sigma | ? | dimensionality contrast | ? |
| n-Transport | A_sigma_dom | ? | dominant-channel anisotropy | ? |
| n-Transport | kappa_e | NO / sensitivity | validation | redundancy |
| n-Transport | PF | NO | external label | derived property |
| p-Transport | ... | ... | ... | ... |
| Electronic | Eg | ? | primary | ? |
| Electronic | m*_median | ? | primary | ? |
| Electronic | m*_spread | ? | spectrum | ? |

所有问号必须根据真实结果填写。

---

# 六十一、最终冻结文件

只有经过审计并确定后，

才能正式写出：

`FINAL_TRANSPORT_REPRESENTATION_V1`

例如最终可能是：

```text
n-Transport V1

S_median
S_MAD
S_sign_fraction
log_sigma_dom_geo
D_sigma
A_sigma_dom
```

以及：

```text
p-Transport V1

S_median
S_MAD
S_sign_fraction
log_sigma_dom_geo
D_sigma
A_sigma_dom
```

这里只是示例。

不要提前固定。

最终必须由 Phase L0 真实结果决定。

---

# 六十二、结果解释的语言规范

没有 eigenvectors 的情况下：

允许：

- principal value；
- principal-value spectrum；
- dominant channel；
- weakest channel；
- suppressed-channel contrast；
- dominant-channel anisotropy；
- quasi-2D-like spectral signature。

禁止直接写：

- x conductivity；
- y conductivity；
- z conductivity；
- in-plane eigenvalue；
- out-of-plane eigenvalue。

除非后续获得 eigenvectors 或原始 tensor。

---

# 六十三、禁止过度解释相关性

例如发现：

Seebeck sign violation 与 small gap 显著相关。

可以写：

> Sign-inconsistent materials are enriched in the small-gap region.

不能直接写：

> Small band gap causes the sign violation.

除非进一步具有因果或物理机制证据。

---

# 六十四、代码质量要求

新增脚本建议：

```text
scripts/
    11_pf_pairing_audit.py
    12_transport_spectrum_audit.py
    13_seebeck_sign_analysis.py
    14_effective_mass_spectrum.py
    15_transport_representation_compare.py
    16_freeze_transport_features.py
    17_freeze_electronic_features.py
```

所有脚本必须：

- 可独立运行；
- 有清楚日志；
- 不静默忽略错误；
- 不覆盖原始数据；
- 对 NaN / Inf 显式检查；
- 输出随机种子；
- 输出使用参数。

---

# 六十五、所有数据必须保留 JID

任何：

CSV

Parquet

JSON

必须保留：

`jid`

作为材料身份锚点。

后续所有不同 View 之间的对齐都依赖 JID。

禁止在中间处理中丢失 JID。

---

# 六十六、本阶段禁止的工作

Phase L0 严格禁止：

- SOAP feature generation；
- Structure kNN graph；
- UMAP；
- t-SNE；
- Diffusion Map；
- Spectral Embedding；
- SNF；
- multilayer graph；
- supra adjacency matrix；
- manifold alignment；
- joint Laplacian；
- candidate superlattice generation。

本阶段只做：

**transport representation physics audit**

以及：

**final feature definition**

---

# 六十七、停止条件

完成以下全部内容后：

- L0-A PF pairing ambiguity；
- L0-B conductivity / κ_e spectrum；
- L0-C Seebeck spectrum；
- L0-D representation comparison；
- L0-E effective-mass spectrum；
- L0-F T1/T2/T3 comparison；
- L0-G n/p Transport feature freeze；
- L0-H Electronic feature freeze；

立即：

**STOP**

不要继续下一阶段。

---

# 六十八、最终返回给我的结果格式

最终回答必须首先给一个 Executive Summary。

然后依次给：

## 1. PF pairing

报告：

- n-type median ambiguity；
- p-type median ambiguity；
- >20% ambiguity fraction；
- PF 是否可以继续使用。

## 2. Conductivity spectrum

报告：

- D_sigma distribution；
- A_sigma_dom distribution；
- quasi-2D-like spectral pattern 是否明显。

## 3. Seebeck

报告：

- N3/N2/N1/N0；
- P3/P2/P1/P0；
- sign consistency 与 Eg 的关系；
- mean vs median 哪一个更稳健。

## 4. Representation comparison

报告：

- R1/R2/R3 distance correlation；
- kNN overlap；
- 推荐版本。

## 5. κ_e redundancy

报告：

- T1/T2/T3 overlap；
- κ_e 是否保留在主模型。

## 6. Effective mass

报告：

- 数据定义；
- spectrum feature；
- n/p 处理方式。

## 7. Final Feature Decision Table

明确列出：

- KEEP；
- DROP；
- SENSITIVITY；
- EXTERNAL LABEL。

## 8. Frozen V1

明确写出最终：

`n_transport_features_v1`

`p_transport_features_v1`

`electronic_features_v1`

具体字段。

完成以上内容后 STOP。

---

# 六十九、下一阶段预告

只有 Phase L0 通过以后，下一阶段才执行：

```text
Phase L
Structure descriptor construction
        ↓
Phase M
Structure similarity graph
        ↓
Phase N
Electronic / n-Transport / p-Transport graphs
        ↓
Phase O
Single-view geometry validation
        ↓
Phase P
Partial multilayer JID alignment
        ↓
Phase Q
Joint spectral / diffusion embedding
        ↓
Phase R
Unified 2D Thermoelectric Transport Atlas
        ↓
Phase S
PF external performance mapping
        ↓
Phase T
Structure-close / Transport-far pair discovery
```

当前绝对不要提前执行这些阶段。

---

# 七十、整个 Phase L0 最重要的科学目标

当前不是为了把特征做得越多越好。

真正目标是：

**从 JARVIS 的 3 个 tensor principal values 中提取最少但物理意义最清晰的输运自由度。**

特别需要区分：

`overall spectrum spread`

和：

`suppressed-channel dimensionality contrast`

以及：

`dominant-channel anisotropy`

不能把三者混成一个简单的：

`max / min`

指标。

同时必须解决：

**Seebeck mean 是否受到异常 principal channel 污染。**

并确认：

**PF 的 eigenvalue pairing convention 是否会影响材料性能排名。**

最终得到的 V1 Transport Representation 必须满足：

- permutation invariant；
- physically interpretable；
- robust to one abnormal principal value；
- preserves meaningful anisotropy；
- avoids σ / κ_e double counting；
- excludes PF target leakage；
- suitable for subsequent manifold construction。

现在开始执行：

**Phase L0：二维输运 principal-value spectrum 物理一致性审计。**

完成以后 STOP，并返回全部结果和最终冻结的 V1 特征集合。