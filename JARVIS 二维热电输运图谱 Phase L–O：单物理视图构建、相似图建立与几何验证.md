# JARVIS 二维热电输运图谱 Phase L–O：单物理视图构建、相似图建立与几何验证

## 一、当前任务定位

前序阶段已经完成：

**Phase A–F**

JARVIS 二维数据库下载、来源验证、Schema Audit 和 Property Coverage Audit。

**Phase G–K**

输运三本征值恢复、JID 交集审计、数值审计和输运变量冗余分析。

**Phase L0**

二维输运 principal-value spectrum 物理一致性审计，并冻结最终 V1 representation。

当前已经明确：

- 不使用 κ_L；
- 不构建 ZT；
- PF 不进入主输运 embedding；
- κ_e 只作为 sensitivity；
- 保留 principal-value spectrum 信息；
- 不把三个 eigenvalues 解释为 x/y/z；
- n-type 和 p-type 分开研究。

当前阶段开始执行：

**Phase L–O：Single-View Geometry Construction and Validation**

本阶段目标是：

> 分别构建 Structure、Electronic、n-Transport、p-Transport 的材料相似图，确认每个物理空间本身具有稳定、合理和可解释的局部几何结构。

本阶段结束以前：

**禁止进行不同 View 的正式融合。**

---

# 二、当前冻结数据

数据库总材料数：

`N = 1103`

---

## 2.1 Structure

覆盖：

`1103 / 1103`

当前还没有正式 structure descriptor。

本阶段建立。

---

## 2.2 n-Transport V1

覆盖：

`806`

冻结特征：

```text
S_median
S_MAD
S_sign_fraction
log_sigma_dom_geo
D_sigma
A_sigma_dom
```

其中：

`S_median`

代表三个 Seebeck principal values 的中位数。

`S_MAD`

描述 Seebeck principal-value spectrum 的稳健离散程度。

`S_sign_fraction`

描述三个 Seebeck 主值中符合该 carrier type 预期符号的比例。

`log_sigma_dom_geo`

描述两个 dominant conductivity principal channels 的几何平均输运尺度。

`D_sigma`

定义为：

`log10(sigma_2 / sigma_1)`

其中：

`sigma_1 ≤ sigma_2 ≤ sigma_3`

表示 weakest principal channel 相对于 dominant channels 的压低程度。

`A_sigma_dom`

定义为：

`log10(sigma_3 / sigma_2)`

描述两个 dominant principal channels 之间的谱各向异性。

---

## 2.3 p-Transport V1

覆盖：

`803`

使用与 n-type 完全相同的 6 个 feature definitions：

```text
S_median
S_MAD
S_sign_fraction
log_sigma_dom_geo
D_sigma
A_sigma_dom
```

但使用 p-type 数据。

---

## 2.4 Electronic V1

Band gap：

`Eg_optb88vdw`

覆盖：

`1103`

Effective mass：

覆盖：

`678`

当前冻结：

```text
m_elec_median
m_hole_median
```

因此 Electronic View 天然存在两种数据层级：

### Electronic-gap layer

覆盖 1103：

```text
Eg_optb88vdw
```

### Electronic-rich layer

覆盖 678：

```text
Eg_optb88vdw
m_elec_median
m_hole_median
```

本阶段禁止通过填补 effective mass，使其人为扩展到 1103。

---

# 三、本阶段核心原则

## 原则 1：先建立每个 View 自己的几何，再考虑融合

当前不要构造：

- supra graph；
- joint Laplacian；
- joint embedding；
- unified materials atlas。

首先确认：

`G_structure`

`G_electronic`

`G_n_transport`

`G_p_transport`

各自是合理的。

---

## 原则 2：不得使用 PF 调参

PF 是后续 external performance label。

禁止为了让 high-PF 材料看起来聚集，而调整：

- SOAP 参数；
- kNN 的 k；
- graph kernel；
- Structure fusion weight；
- transport scaler；
- UMAP 参数。

PF 只能在图冻结之后用于 post-hoc validation。

---

## 原则 3：不得使用最终热电目标构造 Structure View

Structure View 只能使用：

- composition；
- atomic geometry；
- local environment；
- 2D geometry。

禁止输入：

- band gap；
- effective mass；
- Seebeck；
- conductivity；
- PF；
- κ_e。

---

## 原则 4：二维真空层不能成为结构相似性的来源

两个完全相同的二维层：

一个设置：

`15 Å vacuum`

另一个：

`30 Å vacuum`

必须在 Structure View 中几乎完全相同。

---

# 四、本阶段执行流程

严格执行：

```text
Phase L
2D structure preprocessing
        ↓
Structure descriptor construction
        ↓
Structure descriptor invariance tests
        ↓

Phase M
Structure similarity graph
        ↓
Structure graph validation
        ↓

Phase N
Electronic graph
n-Transport graph
p-Transport graph
        ↓

Phase O
Single-view geometry diagnostics
Cross-view neighborhood comparison
External-label diagnostics
        ↓
Freeze all single-view graphs
        ↓
STOP
```

完成 Phase O 后停止。

下一阶段再执行正式：

`Partial Multi-View Alignment`

---

# 五、Phase L：二维结构预处理

## L1. 输入数据

使用：

`data/raw/jarvis/dft_2d_snapshot.json`

或已经验证过的对应 structure table。

所有记录必须保留：

`jid`

作为唯一材料身份。

---

# 六、L2：识别二维非周期方向

不能默认所有结构的第三个 lattice vector 都是 vacuum direction。

需要实际检查。

对于每个材料：

获取：

- lattice；
- fractional coordinates；
- Cartesian coordinates。

对每个 fractional axis：

1. 将坐标映射到 `[0,1)`；
2. 排序；
3. 计算相邻 fractional coordinates 之间的 gap；
4. 包括周期边界：

`last → first + 1`

5. 找到 largest empty fractional gap；
6. 乘以对应 lattice vector length，得到近似 physical empty gap。

将 largest physical empty gap 最大的 axis 作为：

`vacuum-axis candidate`

---

# 七、L3：记录 vacuum-axis confidence

对于三个方向的最大 empty gap：

排序：

`gap_1 ≥ gap_2 ≥ gap_3`

定义：

`vacuum_confidence = gap_1 / (gap_2 + epsilon)`

不要直接用该指标自动删除结构。

输出：

```text
data/audit/vacuum_axis_audit.csv
```

至少包含：

```text
jid
formula
vacuum_axis
vacuum_gap
second_largest_gap
vacuum_confidence
status
```

---

# 八、L4：人工审计异常结构

重点列出：

- vacuum gap 很小；
- vacuum direction 不明显；
- 第一、第二候选方向差异很小；
- 原子跨周期边界导致 slab 看起来被拆开的结构。

输出：

```text
reports/vacuum_axis_anomalies.md
```

如果大量结构无法可靠识别二维方向：

**STOP**

不要继续 SOAP。

如果只是极少数异常：

保留 JID 并单独记录。

---

# 九、L5：统一 2D 周期条件

确定 vacuum axis 后：

将晶格轴重新排列，使 vacuum axis 成为第三方向。

然后对于用于局域结构 descriptor 的 ASE Atoms：

设置：

```python
pbc = [True, True, False]
```

即：

- 两个面内方向周期；
- vacuum direction 非周期。

不要让 SOAP 看到跨真空层的周期镜像。

---

# 十、L6：slab recentering

沿非周期方向将整个 slab 平移到 cell 中心。

注意：

这只是坐标规范化。

不能改变：

- 原子间距离；
- 键；
- 面内结构。

保存：

```text
structures/standardized_2d/
```

同时保留原始结构。

不要覆盖原始数据。

---

# 十一、L7：建立三个 Structure information blocks

Structure View 不要一开始把所有 descriptor 混在一起。

首先分别建立：

### Block G

Geometry。

### Block C

Composition。

### Block S

Species-sensitive local environment。

其中 Block S 为 sensitivity。

主模型首先依赖：

`Geometry + Composition`

再检查 species-sensitive descriptor 是否带来额外信息。

---

# 十二、L8：Composition Block

为了避免引入额外数据库的元素属性和缺失值，

第一版使用最干净的：

**elemental fraction vector**

设数据库全部出现的元素集合为：

`E = {E1, E2, ..., Em}`

对于材料 i：

构造：

```text
c_i =
[fraction(E1),
 fraction(E2),
 ...
 fraction(Em)]
```

所有分量之和：

`sum(c_i) = 1`

---

# 十三、L9：Composition Distance

由于 composition vector 本质上是概率分布，

优先使用 Hellinger distance：

`d_comp(i,j) = sqrt(0.5 × Σ( sqrt(c_i,k) - sqrt(c_j,k) )² )`

优点：

- 有界；
- 对 composition fraction 有明确数学意义；
- 不需要外部元素 property；
- 不存在人为元素编码顺序问题。

同时保留：

cosine distance

作为 sensitivity baseline。

---

# 十四、L10：Geometry SOAP Baseline

第一版 SOAP 首先建立：

**geometry-only SOAP**

方法：

复制结构。

将所有原子临时映射成同一种 dummy species。

例如：

`X`

仅保留：

- atomic positions；
- local geometry；
- coordination environment；
- topology。

这样 SOAP 不会因为数据库元素种类很多而产生非常巨大的 species-channel 维度。

注意：

只是计算 descriptor 时映射。

原始 structure 绝对不能修改。

---

# 十五、L11：SOAP 参数不能随意选择

必须显式记录：

```text
r_cut
n_max
l_max
sigma
periodic setting
averaging method
dtype
```

首先建立少量候选。

例如：

### SOAP-A

较短局域环境：

`r_cut ≈ 4 Å`

### SOAP-B

中等局域环境：

`r_cut ≈ 6 Å`

### SOAP-C

更长局域环境：

`r_cut ≈ 8 Å`

实际值需要结合二维材料典型近邻距离检查。

不要让：

`r_cut`

超过 vacuum gap，从而重新引入跨层镜像。

---

# 十六、L12：Global SOAP representation

不要把所有 atomic SOAP vectors 直接拼接。

对于每个材料：

先计算每个原子的 SOAP vector：

`x_atom`

然后进行 permutation-invariant pooling。

Baseline：

`x_structure = mean(x_atom over atoms)`

可以同时测试：

`mean + std`

但不要直接把每个原子按 POSCAR 顺序展开。

最终每个材料必须对应固定长度的结构 vector。

---

# 十七、L13：SOAP Feature Normalization

对每个 structure SOAP vector：

进行 L2 normalization：

`x_norm = x / ||x||`

然后定义 SOAP kernel similarity：

`K_SOAP(i,j) = dot(x_i, x_j)`

定义 SOAP distance：

`d_SOAP(i,j) = sqrt(max(0, 2 - 2 × K_SOAP(i,j)))`

不要未经验证直接使用原始高维 SOAP 的未归一化 Euclidean distance。

---

# 十八、L14：Species-Sensitive SOAP Sensitivity

如果当前安装的 DScribe 支持可控制维度的 species compression：

检查当前安装版本真实 API。

使用：

```python
inspect.signature(SOAP)
```

并阅读本地 documentation。

不要凭记忆猜参数名。

如果存在可靠的 compressed species-sensitive SOAP：

建立：

`SOAP_species`

作为 sensitivity descriptor。

如果 species-sensitive SOAP 维度极大：

例如：

- 占用内存不可接受；
- feature dimension 远大于合理范围；

不要硬算。

记录：

`SKIPPED_HIGH_DIMENSION`

Baseline 仍采用：

`geometry-only SOAP + composition`

---

# 十九、L15：二维 Structure Invariance Test

这是 SOAP 正式进入数据库分析前的硬性 QA。

随机选择至少：

`50 个材料`

进行以下测试。

---

# 二十、Test A：Atom permutation invariance

随机改变原子排列顺序。

结构不变。

重新计算 descriptor。

要求：

SOAP similarity 极接近：

`1`

working tolerance：

`> 0.999999`

---

# 二十一、Test B：Slab translation invariance

将整个二维层沿：

- 面内方向；
- 非周期方向；

整体平移。

不改变内部结构。

重新计算 descriptor。

SOAP similarity 应接近：

`1`

---

# 二十二、Test C：Vacuum invariance

对于同一材料，

构造：

`15 Å`

`20 Å`

`25 Å`

`30 Å`

真空版本。

原子层结构不变。

重新计算 geometry descriptor。

要求不同 vacuum 版本之间：

SOAP similarity：

`> 0.999`

working threshold。

如果明显低于该值：

当前结构 representation 不合格。

停止并修正。

---

# 二十三、Test D：Supercell invariance

随机选择至少：

`20`

个材料。

构造面内：

`2 × 2`

supercell。

使用 global mean SOAP pooling。

比较 primitive/canonical representation 与 supercell representation。

理论上局域环境统计应该高度一致。

目标：

SOAP similarity：

`> 0.995`

如果明显偏低：

检查：

- cutoff；
- pooling；
- PBC；
- cell preprocessing。

---

# 二十四、L16：Structure QA 输出

生成：

```text
data/audit/structure_descriptor_invariance.csv
reports/structure_descriptor_invariance.md
```

只有所有主要 invariance tests 通过：

才能进入 Phase M。

---

# 二十五、Phase M：Structure Similarity Graph

首先建立三个候选 Structure similarities。

---

# 二十六、M1：Geometry Graph

基于：

`d_SOAP`

构建：

`G_geo`

---

# 二十七、M2：Composition Graph

基于：

`d_comp`

构建：

`G_comp`

---

# 二十八、M3：Combined Structure Graph

不要把 SOAP vector 与 composition vector 直接粗暴拼接。

首先分别形成：

`W_geo`

和：

`W_comp`

再在 affinity level 融合。

Baseline：

`W_struct = 0.5 × W_geo + 0.5 × W_comp`

这是预注册 baseline。

不要根据 PF 调整该比例。

---

# 二十九、M4：Structure Fusion Sensitivity

额外测试：

```text
Geometry : Composition

0.25 : 0.75
0.50 : 0.50
0.75 : 0.25
```

如果有 species-sensitive SOAP：

再建立一个 sensitivity structure graph。

目标不是选择让 PF 最漂亮的结果。

而是检查：

**Structure neighborhood 是否对合理的 chemistry/geometry weighting 稳定。**

---

# 三十、M5：kNN Graph Construction

对于每个 distance matrix：

测试：

```text
k = 10
15
20
30
40
50
```

不要默认 k=10。

---

# 三十一、M6：Local-Scale Affinity

对于材料 i：

定义：

`sigma_i = distance to its k-th nearest neighbor`

对于 kNN edge：

定义：

`W_ij = exp[ -d_ij² / (sigma_i × sigma_j + epsilon) ]`

其中 epsilon 仅用于数值稳定。

应根据浮点数尺度选择非常小的数，并记录。

---

# 三十二、M7：Graph Symmetrization

分别测试：

### Union-kNN

如果：

`i → j`

或：

`j → i`

存在邻居关系，

保留该 edge。

### Mutual-kNN

只有：

`i → j`

且：

`j → i`

同时成立才保留。

Mutual-kNN 更严格，但可能造成 graph fragmentation。

---

# 三十三、M8：选择 k 的规则

禁止根据 property performance 选择 k。

根据纯图结构决定。

至少检查：

- number of connected components；
- giant-component fraction；
- isolated nodes；
- mean degree；
- median degree；
- degree distribution；
- neighbor stability。

优先选择：

**保证绝大多数节点位于同一个 giant component 的最小合理 k。**

working target：

`giant component > 99%`

但这不是不可改变的自然定律。

必须同时报告实际结果。

---

# 三十四、M9：Structure Neighbor Stability

比较：

`k = 10 vs 15`

`15 vs 20`

`20 vs 30`

等。

使用：

Jaccard neighbor overlap。

同时比较：

`0.25/0.75`

`0.5/0.5`

`0.75/0.25`

三种 Structure block weight。

如果结构邻域随微小参数改变而完全重排：

说明 representation 不稳定。

不要继续 joint manifold。

---

# 三十五、M10：Structure Near-Duplicate Audit

寻找：

`d_struct`

非常小的材料对。

例如：

最低：

`0.1%`

距离。

输出前：

`100`

对。

检查：

- same formula？
- same stoichiometry？
- same prototype？
- different JID？
- possible duplicate？

输出：

```text
data/audit/structure_near_duplicates.csv
```

如果 geometry-only SOAP 把完全不同元素但同构结构视为非常近：

这是预期现象。

Combined Structure Graph 应通过 composition block 将其合理区分。

---

# 三十六、M11：冻结 Structure Graph

只有完成上述测试后：

冻结：

`G_structure_v1`

保存：

```text
graphs/G_structure_v1.npz
graphs/G_structure_v1_nodes.csv
graphs/G_structure_v1_metadata.json
```

metadata 必须包含：

- descriptor version；
- SOAP parameters；
- composition distance；
- fusion weights；
- k；
- graph type；
- normalization；
- creation date。

---

# 三十七、Phase N：Electronic Similarity Layers

Electronic View 由于 effective mass 不完整：

不要强行只有一张图。

本阶段建立三个 Electronic candidates。

---

# 三十八、N1：Band-Gap Layer

覆盖：

`1103`

只使用：

`Eg_optb88vdw`

建立：

`G_Eg`

这严格来说是：

**band-gap similarity layer**

不是复杂高维 manifold。

距离：

`|Eg_i - Eg_j|`

经过适当 robust scale normalization。

不要把金属：

`Eg = 0`

当 missing。

---

# 三十九、N2：n-Electronic Rich Layer

覆盖具有 effective mass 的：

`678`

使用：

```text
Eg_optb88vdw
m_elec_median
```

建立：

`G_electronic_n`

其目标是与 n-type transport 保持 carrier physics 对应。

---

# 四十、N3：p-Electronic Rich Layer

同样在 678 个材料上使用：

```text
Eg_optb88vdw
m_hole_median
```

建立：

`G_electronic_p`

---

# 四十一、N4：Generic Electronic Layer

作为 sensitivity：

```text
Eg_optb88vdw
m_elec_median
m_hole_median
```

覆盖：

`678`

建立：

`G_electronic_joint`

这不是默认主层。

它主要用于判断：

将 electron/hole mass 混在一起是否改变电子材料邻域。

---

# 四十二、N5：Electronic Scaling

必须读取：

`features/electronic/electronic_feature_metadata.csv`

遵循 Phase L0 已冻结的 property definitions。

在进入距离计算前：

检查：

- distribution；
- outliers；
- positive/negative；
- heavy tail。

推荐 baseline：

`RobustScaler`

但不能改变 feature 的物理定义。

如果 m* 严重跨数量级：

可以测试 log transform。

但必须：

1. 只对物理上允许且严格正的数据；
2. 记录 transformation；
3. 与 raw/robust baseline 比较；
4. 不根据 PF 决定使用哪一个。

---

# 四十三、Phase N：n-Transport Graph

输入：

`features/transport/n_transport_features_v1.parquet`

必须严格使用 Phase L0 冻结的 6 个变量。

不得偷偷重新加入：

- S_mean；
- PF；
- κ_e。

---

# 四十四、N6：n-Transport Scaling

推荐首先测试：

`RobustScaler`

原因：

部分谱离散程度和 dimensionality contrast 可能存在 heavy tail。

同时建立：

`StandardScaler`

sensitivity。

比较二者邻域稳定性。

不要根据 PF 选择 scaler。

---

# 四十五、N7：n-Transport Distance

Baseline：

标准化后的 6 维 feature space：

Euclidean distance。

同时测试：

cosine distance

作为 sensitivity。

如果两者材料邻域高度一致：

优先保留更简单的 Euclidean baseline。

---

# 四十六、N8：Feature Dominance Audit

对 6 个 feature：

分别进行 leave-one-feature-out。

例如：

删除：

`S_median`

重新计算近邻。

删除：

`D_sigma`

重新计算近邻。

……

计算：

`kNN overlap(full, leave-one-out)`

目的：

识别是否某一个 feature 几乎完全控制 transport geometry。

这不是要自动删除 feature。

只是理解每一个 transport degree of freedom 的作用。

---

# 四十七、N9：冻结 n-Transport Graph

使用和 Structure Graph 相同的：

- k scanning；
- connectivity；
- local-scale affinity；
- union/mutual sensitivity；
- graph stability；

流程。

最终冻结：

```text
graphs/G_n_transport_v1.npz
graphs/G_n_transport_v1_nodes.csv
graphs/G_n_transport_v1_metadata.json
```

覆盖：

`806`

---

# 四十八、Phase N：p-Transport Graph

完全使用相同程序处理：

`p_transport_features_v1.parquet`

冻结：

```text
graphs/G_p_transport_v1.npz
graphs/G_p_transport_v1_nodes.csv
graphs/G_p_transport_v1_metadata.json
```

覆盖：

`803`

n/p 两套流程必须尽可能对称。

---

# 四十九、N10：κ_e Sensitivity Graph

根据 Phase L0 的 T2 representation：

分别建立：

`G_n_transport_kappa_sensitivity`

`G_p_transport_kappa_sensitivity`

只用于 sensitivity。

不作为主 View。

比较：

`σ-based graph`

与：

`κ_e-based graph`

的：

- distance-rank correlation；
- kNN overlap；
- graph connectivity。

如果 overlap 很高：

进一步确认：

κ_e 不需要进入主模型。

---

# 五十、Phase O：Single-View Geometry Validation

Phase O 的核心问题是：

> 每个物理空间本身到底有没有稳定的材料邻域？

不要急着把它们融合。

---

# 五十一、O1：每张图的 Graph QA

对：

`G_structure_v1`

`G_Eg`

`G_electronic_n`

`G_electronic_p`

`G_n_transport_v1`

`G_p_transport_v1`

分别统计：

```text
N nodes
N edges
number of connected components
giant component fraction
mean degree
median degree
min degree
max degree
isolated nodes
edge-weight distribution
```

保存：

```text
data/audit/single_view_graph_QA.csv
```

---

# 五十二、O2：Single-View Spectral Diagnostic

现在才允许第一次进行 spectral decomposition。

对于每一个冻结 graph：

构造 normalized adjacency 或 normalized graph Laplacian。

例如：

`D_ii = Σ_j W_ij`

然后：

`L_sym = I - D^(-1/2) W D^(-1/2)`

求最小的一组 non-trivial eigenvalues/eigenvectors。

主要目的是：

- 查看 eigenvalue spectrum；
- 查看 eigengap；
- 初步判断 manifold complexity；
- 检查 graph 是否存在人为断裂。

---

# 五十三、O3：Single-View Embedding

允许计算：

### Spectral Embedding

作为主要几何诊断。

### Diffusion-style coordinates

如果已经有稳定实现，可以计算。

### UMAP

仅作为 visual diagnostic。

UMAP 不能代替 graph distance。

---

# 五十四、O4：UMAP 不参与 Graph 选择

禁止：

> 哪个 k 的 UMAP 图看起来最漂亮，就选哪个 k。

k 和 graph parameters 必须在 UMAP 之前依据：

- graph connectivity；
- neighbor stability；
- representation robustness；

冻结。

UMAP 只能在冻结以后画。

---

# 五十五、O5：Cross-View Neighbor Overlap

虽然当前不做 alignment，

现在可以比较不同 View 原始邻域是否一致。

对于两个 View：

只取共同存在的 JIDs。

例如：

`Structure ∩ n-Transport`

计算：

`NN_structure(i)`

和：

`NN_n_transport(i)`

然后：

`Overlap_i = |NN_A(i) ∩ NN_B(i)| / k`

最终平均。

---

# 五十六、O6：必须比较的 View Pairs

至少计算：

```text
Structure vs Eg

Structure vs n-Electronic

Structure vs p-Electronic

Structure vs n-Transport

Structure vs p-Transport

n-Electronic vs n-Transport

p-Electronic vs p-Transport

n-Transport vs p-Transport
```

其中最后一个只在：

共同 n/p JIDs

上计算。

---

# 五十七、O7：Cross-View Distance Correlation

对于共同 JIDs：

计算两个 View 的 pairwise material distances。

然后计算：

`Spearman(distance_A, distance_B)`

为了节省内存：

如果 pair 数过大，可以固定随机种子采样：

最多：

`200000`

个 material pairs。

必须记录采样方法和 random seed。

---

# 五十八、O8：Neighbor-Overlap Random Baseline

不能只看到：

`neighbor overlap = 0.18`

就判断高或低。

需要建立随机 baseline。

随机打乱一个 View 的 JID mapping：

`1000 次`

重新计算 neighbor overlap。

得到：

- null mean；
- null std；
- percentile；
- empirical p-value。

这样才能判断：

不同物理空间的邻域重合是否显著高于随机。

---

# 五十九、O9：建立 View Similarity Matrix

最终得到：

```text
                 Struct   Eg   E-n   E-p   T-n   T-p

Structure

Eg

Electronic-n

Electronic-p

Transport-n

Transport-p
```

矩阵元素至少包括：

### Matrix A

mean kNN overlap。

### Matrix B

distance-rank correlation。

输出：

```text
data/audit/view_neighbor_overlap.csv
data/audit/view_distance_correlation.csv
```

图：

```text
figures/view_neighbor_overlap_heatmap.png
figures/view_distance_correlation_heatmap.png
```

---

# 六十、O10：PF External Label Diagnostics

到这里所有主 graph 已经冻结。

现在才允许使用 PF。

PF 不参与 graph construction。

---

# 六十一、O11：分别处理 n-PF 和 p-PF

使用：

`PF_n`

只分析：

`G_n_transport_v1`

以及共同 JID 的：

- Structure；
- Electronic-n。

使用：

`PF_p`

只分析：

`G_p_transport_v1`

以及：

- Structure；
- Electronic-p。

n/p 禁止混合。

---

# 六十二、O12：PF 数值定义说明

Phase L0 已确认：

PF 对 eigenvalue pairing 高度敏感。

因此所有报告必须使用：

**JARVIS-convention database-defined PF**

或者简称：

**JARVIS-defined PF**

禁止写成：

**true directional PF**

也不要声称它对应某一个 crystallographic direction。

---

# 六十三、O13：PF Graph Smoothness

对于 frozen graph W 和 PF：

计算：

`Smoothness(PF) = Σ W_ij × (PF_i - PF_j)² / Σ W_ij`

建议对 PF 先考虑：

`log10(PF + epsilon)`

因为 PF 往往跨多个数量级。

epsilon 需要根据正值数据确定并记录。

---

# 六十四、O14：PF Smoothness Null Test

随机打乱 PF：

`1000 次`

每次重新计算 smoothness。

比较真实 smoothness：

与随机 null distribution。

输出：

- empirical p-value；
- z-score；
- percentile。

这回答：

> 高低 PF 是否在 transport graph 上具有局域连续性？

注意：

这是 post-hoc validation。

不是 graph tuning。

---

# 六十五、O15：不同 View 的 PF Smoothness

分别比较：

n-type：

```text
Structure graph
Electronic-n graph
n-Transport graph
```

上的 PF smoothness。

p-type 同理。

如果：

PF 在 Transport Graph 上明显更平滑，

说明 transport representation 捕获了与 PF 有关的基础物理信息。

如果 PF 在 Structure Graph 上也很平滑：

说明结构本身可能已经强烈约束输运。

---

# 六十六、O16：Band-Gap External Mapping

在 Structure Graph 冻结以后：

可以使用 Eg 作为外部验证量。

检查：

`Eg`

在 Structure Graph 上是否平滑。

目的：

不是证明 Structure Graph 好坏的唯一标准，

而是判断：

> 结构局部邻域是否与电子结构存在统计连续关系。

---

# 六十七、O17：n/p Transport Relationship

对同时具有 n 和 p transport 的材料：

比较：

`G_n_transport`

与：

`G_p_transport`

回答：

> 一个材料的 n-type transport 邻域和 p-type transport 邻域有多一致？

计算：

- neighbor overlap；
- distance-rank correlation。

如果很低：

说明 n/p transport 是两个真正不同的材料空间。

这会支持后续将两者作为独立 View。

---

# 六十八、O18：Structure–Transport Tension Preview

当前还没有 joint embedding，

因此不要使用之前定义的 joint-space cross-view tension。

但可以定义一个 preliminary neighborhood disagreement：

`Disagreement_i = 1 - Overlap_i(structure, transport)`

分别：

`Disagreement_n`

`Disagreement_p`

---

# 六十九、O19：寻找高 disagreement 材料

列出：

Structure vs n-Transport：

disagreement 最高前 50 个材料。

Structure vs p-Transport：

同样前 50 个。

保存：

```text
data/audit/high_structure_transport_disagreement_n.csv
data/audit/high_structure_transport_disagreement_p.csv
```

这里只作为 diagnostic。

还不能正式称为：

`structure-close / transport-far superlattice candidates`

那属于后续 Phase T。

---

# 七十、O20：人工检查 Neighbor Examples

对于每个主 View：

随机抽取至少：

`20`

个 anchor materials。

列出其：

`10 nearest neighbors`

输出：

```text
reports/neighbor_examples_structure.md
reports/neighbor_examples_electronic_n.md
reports/neighbor_examples_electronic_p.md
reports/neighbor_examples_transport_n.md
reports/neighbor_examples_transport_p.md
```

每条显示：

- JID；
- formula；
- distance；
- relevant feature values。

人工检查是否出现明显荒谬邻居。

---

# 七十一、Structure Neighbor 特别检查

例如：

两个材料：

composition 完全不同，

但 geometry 相同。

如果 Geometry-only SOAP 把它们放得很近：

正常。

但 Combined Structure Graph 应该根据 composition 对它们有所区分。

反之：

同 composition + 相似 structure

却距离极远，

需要检查 descriptor。

---

# 七十二、Transport Neighbor 特别检查

两个 n-type transport 近邻应在 V1 feature 上表现出整体相似：

- S_median；
- spectrum spread；
- conductivity scale；
- dimensionality contrast；
- dominant-channel anisotropy。

不能只因为单个 feature 数值接近就成为近邻。

---

# 七十三、O21：Feature Leave-One-Out Stability

对 n/p Transport View 的 6 个 frozen features：

记录每个 feature 被删除后：

`median kNN overlap`

建立：

```text
data/audit/transport_feature_leave_one_out.csv
```

如果某 feature 删除以后：

邻域几乎完全变化，

说明这个 feature 是主要 geometry driver。

报告。

不要自动删除。

---

# 七十四、O22：Structure Block Ablation

分别比较：

```text
Geometry only

Composition only

Geometry + Composition
```

三种 Structure Graph。

计算：

- neighbor overlap；
- distance-rank correlation；
- graph connectivity。

这将回答：

> 当前二维材料图主要由 chemistry 驱动，还是由 geometry 驱动？

这是后续论文非常有价值的结果。

---

# 七十五、O23：不得使用 PF 选择 Structure 权重

即使发现：

`Geometry 0.75 + Composition 0.25`

对应 PF smoothness 更高，

也不能因为这一点选择该权重。

Structure fusion weight 必须根据：

- graph stability；
- chemistry/geometry 平衡；
- unsupervised robustness；

确定。

PF 结果只能事后报告。

---

# 七十六、O24：最终冻结 Single-View Graphs

完成所有 QA 后，

冻结以下 graph versions。

至少：

```text
G_structure_v1

G_Eg_v1

G_electronic_n_v1

G_electronic_p_v1

G_n_transport_v1

G_p_transport_v1
```

以及 sensitivity：

```text
G_electronic_joint_sensitivity

G_n_transport_kappa_sensitivity

G_p_transport_kappa_sensitivity
```

---

# 七十七、所有 Graph Metadata 必须完整

每张 graph 对应一个 JSON：

```text
graph_name
view
N_nodes
features
feature_transform
distance_metric
k
kernel
kernel_parameters
symmetrization
normalization
connected_components
giant_component_fraction
random_seed
creation_date
```

以后任何 joint analysis 必须读取 metadata，

不能凭记忆重建 graph。

---

# 七十八、本阶段建议新增脚本

```text
scripts/
    18_audit_2d_geometry.py
    19_standardize_2d_structures.py
    20_build_structure_descriptors.py
    21_test_structure_invariance.py
    22_build_structure_graph.py
    23_build_electronic_graphs.py
    24_build_transport_graphs.py
    25_single_view_graph_qa.py
    26_single_view_spectral_diagnostics.py
    27_cross_view_neighbor_analysis.py
    28_external_label_diagnostics.py
    29_freeze_single_view_graphs.py
    30_phase_LO_summary.py
```

---

# 七十九、本阶段核心文件输出

## Structure

```text
features/structure/composition_fraction.parquet
features/structure/geometry_soap_v1.parquet
features/structure/structure_feature_metadata.json
```

如有 species-sensitive SOAP：

```text
features/structure/species_soap_sensitivity.parquet
```

---

## Graphs

```text
graphs/G_structure_v1.npz
graphs/G_Eg_v1.npz
graphs/G_electronic_n_v1.npz
graphs/G_electronic_p_v1.npz
graphs/G_n_transport_v1.npz
graphs/G_p_transport_v1.npz
```

---

## Audit

```text
data/audit/vacuum_axis_audit.csv
data/audit/structure_descriptor_invariance.csv
data/audit/structure_near_duplicates.csv
data/audit/single_view_graph_QA.csv
data/audit/view_neighbor_overlap.csv
data/audit/view_distance_correlation.csv
data/audit/transport_feature_leave_one_out.csv
data/audit/high_structure_transport_disagreement_n.csv
data/audit/high_structure_transport_disagreement_p.csv
```

---

# 八十、图形输出

至少生成：

1. Structure graph degree distribution。

2. n-Transport graph degree distribution。

3. p-Transport graph degree distribution。

4. Structure spectral eigenvalue spectrum。

5. n-Transport spectral eigenvalue spectrum。

6. p-Transport spectral eigenvalue spectrum。

7. Structure single-view spectral map。

8. n-Transport single-view spectral map。

9. p-Transport single-view spectral map。

10. View neighbor-overlap heatmap。

11. View distance-correlation heatmap。

12. PF smoothness comparison。

13. Structure block ablation comparison。

14. n/p transport neighbor comparison。

所有图：

- 白底；
- 单图文件；
- 标签完整；
- 保存 PNG；
- 同时保存 PDF 或 SVG 矢量版本。

---

# 八十一、本阶段报告

生成：

```text
reports/structure_preprocessing.md
reports/structure_descriptor_validation.md
reports/structure_graph_selection.md
reports/electronic_graph_analysis.md
reports/transport_graph_analysis.md
reports/cross_view_geometry_analysis.md
reports/external_label_diagnostics.md
reports/phase_LO_summary.md
```

---

# 八十二、Phase O 最终必须回答的问题

完成以后，必须明确回答以下问题。

1. 1103 个二维结构中，多少能够明确识别 vacuum axis？

2. 是否存在异常二维结构？

3. SOAP 是否通过：
   - atom permutation；
   - translation；
   - vacuum；
   - supercell；
   四类 invariance tests？

4. Geometry-only 和 Composition-only 的材料邻域有多大差异？

5. Combined Structure Graph 是否稳定？

6. 最终 Structure Graph 使用什么 SOAP 参数？

7. 最终 Structure Graph 使用什么 geometry/composition weight？

8. 最终 k 是多少？

9. 为什么选择这个 k？

10. Structure Graph giant component fraction 是多少？

11. n-Transport Graph 是否连通和稳定？

12. p-Transport Graph 是否连通和稳定？

13. 哪些 transport features 对邻域影响最大？

14. κ_e sensitivity graph 与 σ-based 主图有多相似？

15. Eg-only layer 与 rich Electronic layer 差别多大？

16. Electron-mass Electronic View 与 n-Transport 有多高 neighbor overlap？

17. Hole-mass Electronic View 与 p-Transport 有多高 neighbor overlap？

18. Structure 与 n-Transport 的 neighbor overlap 是多少？

19. Structure 与 p-Transport 的 neighbor overlap 是多少？

20. n-Transport 与 p-Transport 是否是明显不同的材料空间？

21. 所有真实 cross-view overlaps 是否显著高于随机 baseline？

22. PF 在 n-Transport Graph 上是否具有显著 smoothness？

23. PF 在 p-Transport Graph 上是否具有显著 smoothness？

24. PF 在 Structure Graph 和 Transport Graph 上哪一个更连续？

25. 哪些材料具有最大的 Structure–Transport neighborhood disagreement？

---

# 八十三、最终 Single-View Decision Table

生成：

| Layer | N | Features | Distance | k | Status | Role |
|---|---:|---|---|---:|---|---|
| Structure | 1103 | SOAP + composition | ? | ? | ? | backbone |
| Eg | 1103 | Eg | ? | ? | ? | electronic scalar layer |
| Electronic-n | 678 | Eg + m_elec | ? | ? | ? | n electronic |
| Electronic-p | 678 | Eg + m_hole | ? | ? | ? | p electronic |
| n-Transport | 806 | frozen V1 | ? | ? | ? | primary |
| p-Transport | 803 | frozen V1 | ? | ? | ? | primary |
| κe-n | 806 | sensitivity | ? | ? | sensitivity | secondary |
| κe-p | 803 | sensitivity | ? | ? | sensitivity | secondary |

所有问号必须根据真实运行结果填写。

---

# 八十四、必须明确区分 Conceptual View 与 Computational Layer

论文概念上仍然可以描述为四大物理视图：

```text
Structure

Electronic

n-Transport

p-Transport
```

但是计算实现中：

Electronic 可以包含：

```text
Eg layer
Electronic-n rich layer
Electronic-p rich layer
```

因为 effective mass 数据不完整。

因此不要为了维持“四张图”的表面整洁而强行补全数据。

---

# 八十五、本阶段禁止的工作

严格禁止：

- supra adjacency matrix；
- JID anchor coupling；
- multilayer graph alignment；
- joint spectral embedding；
- joint Diffusion Map；
- unified materials atlas；
- high-PF ridge detection；
- superlattice candidate ranking；
- ML prediction；
- missing-value imputation。

---

# 八十六、停止条件

只有当：

`G_structure_v1`

`G_electronic_n_v1`

`G_electronic_p_v1`

`G_n_transport_v1`

`G_p_transport_v1`

以及必要辅助 graph 全部完成：

- QA；
- stability；
- connectivity；
- neighbor analysis；
- external-label diagnostics；

以后：

**STOP**

不要进入 multi-view alignment。

---

# 八十七、最终返回格式

最终回答首先给：

# Executive Summary

至少总结五项：

### 1. Structure representation

最终使用什么 descriptor、参数和权重。

### 2. Graph parameters

Structure / Electronic / Transport 最终分别使用什么 k。

### 3. Single-view stability

哪些图最稳定，哪些存在问题。

### 4. Cross-view relationship

Structure、Electronic、n/p Transport 之间的 neighbor overlap 和 distance correlation。

### 5. PF external validation

PF 在哪个物理空间最平滑。

---

然后给：

# Final Frozen Graphs

明确列出：

```text
G_structure_v1
G_Electronic_n_v1
G_Electronic_p_v1
G_n_transport_v1
G_p_transport_v1
```

以及各自：

- N；
- k；
- features；
- distance；
- graph normalization。

最后给：

# Ready for Phase P?

明确回答：

`YES`

或者：

`NO`

如果 NO：

说明哪一张单视图 graph 仍有问题。

完成以后 STOP。

---

# 八十八、下一阶段预告

只有本阶段通过以后，

下一阶段才开始真正研究：

```text
Structure Graph
       │
       ├──────── Electronic-n
       │
       ├──────── Electronic-p
       │
       ├──────── n-Transport
       │
       └──────── p-Transport
                 │
                 ▼
        JID identity anchors
                 │
                 ▼
        Partial multilayer graph
                 │
                 ▼
        Joint spectral embedding
                 │
                 ▼
   Unified Thermoelectric Transport Atlas
```

下一阶段将重点解决：

- 各 layer 权重如何定义；
- JID anchor 强度 λ 如何选择；
- 缺失 View 如何自然保留；
- 是否应该让 n/p Electronic 与 n/p Transport 分支对齐；
- Structure backbone 是否会过度主导；
- 如何做 random-anchor negative control；
- 如何定义 joint-space cross-view tension；
- 如何寻找 high-PF region；
- 如何寻找 structure-close / transport-far 材料对。

但当前不要提前执行。

---

# 八十九、本阶段最终科学目标

Phase L0 已经解决：

**每个材料在一个 Transport View 中应该如何表示。**

Phase L–O 要解决的是：

> **这些 representation 是否真的形成稳定、有意义、可重复的材料局部几何？**

只有当单个物理 View 本身成立，

才有资格进一步讨论：

> 不同物理流形之间怎样对齐。

因此当前优先级是：

**Single-view validity > beautiful visualization > joint manifold。**

现在开始执行：

**Phase L：二维 Structure representation 与 invariance audit。**

按照 Phase L → M → N → O 顺序执行。

完成 Phase O 后 STOP。