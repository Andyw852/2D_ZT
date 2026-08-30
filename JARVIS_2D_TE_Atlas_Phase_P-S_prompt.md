# JARVIS 二维热电输运图谱 Phase P–S
## 不完备多视图对齐、公共潜在空间验证与统一 Transport Atlas

## 1. 你的角色

你现在是我的二维热电材料、材料信息学、多视图学习、谱图理论、图机器学习和流形学习研究助手。

我要继续此前已经完成的 JARVIS 二维热电输运图谱项目。

当前阶段的目标不是继续做单一物性图，也不是立即寻找超晶格，而是回答一个更基础的问题：

> Structure、Electronic、n-Transport 和 p-Transport 这些局部几何明显不同的材料空间，是否真的能够通过共同材料身份 JID 对齐到一个稳定的公共潜在空间？

这个问题必须由数据决定。

禁止预设一定存在统一流形。

如果数据不支持公共流形，最终结论可以是 Linked Multi-View Atlas，而不是 Unified Multi-View Atlas。这两种结果都是允许的。

---

# 2. 已完成阶段

## Phase A–F

已完成：
- JARVIS `dft_2d` 数据获取；
- Schema Audit；
- Property Coverage Audit；
- 数据来源验证；
- 数据快照；
- JID 身份体系。

数据库总材料数：`1103`

## Phase G–K

已完成：
- JARVIS 官方 NIST XML 输运数据恢复；
- Seebeck / conductivity / κ_e 三本征值恢复；
- effective mass 三本征值恢复；
- JID overlap；
- 数据数值审计；
- σ 与 κ_e 冗余分析。

## Phase L0

已冻结：

### n-Transport V1
覆盖：`806`

```text
S_median
S_MAD
S_sign_fraction
log_sigma_dom_geo
D_sigma
A_sigma_dom
```

### p-Transport V1
覆盖：`803`

使用相同 6 个特征。

### Electronic V1

```text
Eg_optb88vdw
m_elec_median
m_hole_median
```

Eg 覆盖 1103；effective mass 覆盖 678。

## Phase L–O

已冻结：

```text
G_structure_v1
G_Eg_v1
G_electronic_n_v1
G_electronic_p_v1
G_n_transport_v1
G_p_transport_v1
```

主图统一使用 `k = 15`。

---

# 3. 已确认的重要科学事实

## 3.1 Structure View

Structure representation：

```text
geometry-only SOAP
+
elemental fraction
```

SOAP：

```text
dummy species
n_max = 6
l_max = 6
sigma = 1.0
r_cut = 6 Å
mean pooling
L2 normalization
147 dimensions
```

Composition：81 元素 elemental-fraction vector。

距离融合：

```text
d_structure
=
0.5 × d_geometry
+
0.5 × d_composition
```

Structure graph：
- N = 1103；
- k = 15；
- giant component = 1.0；
- 0 isolated。

## 3.2 Structure 与 property space 几乎去相关

Structure vs Electronic / Transport：

```text
distance correlation ≈ 0.02–0.07
neighbor overlap ≈ 1.5%
```

基本接近随机 baseline。

因此本阶段禁止默认 `Structure = absolute backbone`。

Structure 是完整 View，但不代表它天然应该主导公共潜在空间。

## 3.3 Electronic 与 Transport 有一定公共结构

```text
distance correlation ≈ 0.32–0.38
```

属于中等相关。

## 3.4 n-Transport 与 p-Transport 不等价

```text
neighbor overlap ≈ 10.7%
distance correlation ≈ 0.42
z ≈ 60 relative to random
```

因此 n/p 必须先分别建立 Atlas。

## 3.5 PF 的角色

PF 只能定义为 `JARVIS-defined PF` 或 `JARVIS-convention PF`。

PF：
- 不进入 embedding；
- 不参与 λ 选择；
- 不参与 View 权重选择；
- 只作为 external performance label。

## 3.6 κ_e 的角色

n-type：
```text
Pearson ≈ 0.960
Spearman ≈ 0.928
```

p-type：
```text
Pearson ≈ 0.965
Spearman ≈ 0.947
```

主 Transport View 使用 σ，κ_e 只用于 sensitivity。

---

# 4. 本阶段总体目标

```text
Phase P
Partial Multi-View Graph Construction
        ↓
Phase Q
Joint Manifold Feasibility and Validation
        ↓
Phase R
Frozen Unified / Linked Atlas
        ↓
Phase S
PF External Mapping and n/p Comparison
        ↓
STOP
```

本阶段结束后不得进入超晶格设计或最终 parent-pair ranking。

---

# 5. 核心科学问题

必须回答：

1. n-type 是否存在可接受的公共 latent space？
2. p-type 是否存在？
3. Structure 是否应该进入公共空间？
4. Electronic + Transport 是否比 Structure + Electronic + Transport 更容易形成稳定公共几何？
5. JID identity 是否真正改善跨 View 对齐？
6. 是否存在稳定的 inter-view coupling λ？
7. 对齐后各 View 原始局部邻域能保留多少？
8. 是否存在 alignment–preservation trade-off？
9. PF 是否在最终 joint space 中表现出新的连续结构？
10. 若统一流形失败，Linked Multi-View Atlas 是否更符合真实数据？

---

# 6. 使用 Consensus Identity Node

不再使用 Structure-centered star anchor。

对每一个 JID 创建公共身份节点：

`C_i`

它：
- 不携带 descriptor；
- 不携带 property；
- 只表示材料身份。

如果材料 i 具有 Structure、Eg、electron mass、n-Transport，则创建：

```text
S_i
Eg_i
Me_i
Tn_i
```

并通过 `C_i` 连接：

```text
             S_i
              |
              |
Eg_i ------- C_i ------- Tn_i
              |
              |
             Me_i
```

跨 View 只允许 same-JID identity edge。

---

# 7. Missing View 原则

如果没有某个 property：
- 不创建对应节点；
- 不填 0；
- 不 impute；
- 不 ML 预测。

```text
missing property = missing node
```

---

# 8. Identity Edge 完整度归一化

设材料 i 实际拥有 `m_i` 个 View。

主模型使用：

```text
w_anchor(i,v) = λ / m_i
```

因此：

```text
sum_v w_anchor(i,v) = λ
```

避免高完整度材料被过度加权。

---

# 9. 主模型优先建立两个 Atlas

## n-Type Atlas

```text
Structure
Eg
electron effective mass
n-Transport
```

## p-Type Atlas

```text
Structure
Eg
hole effective mass
p-Transport
```

---

# 10. Phase P-A：建立 mass-only graphs

由于 `G_electronic_n_v1` 和 `G_electronic_p_v1` 已经包含 Eg，为避免 Eg 重复计权，主模型另建：

```text
G_m_electron_v1
G_m_hole_v1
```

仅使用：
- `m_elec_median`
- `m_hole_median`

参数：

```text
RobustScaler
Euclidean distance
k = 15
local-scale affinity
```

输出：

```text
graphs/G_m_electron_v1.npz
graphs/G_m_electron_v1_nodes.csv
graphs/G_m_electron_v1_metadata.json

graphs/G_m_hole_v1.npz
graphs/G_m_hole_v1_nodes.csv
graphs/G_m_hole_v1_metadata.json
```

原 `G_electronic_n_v1` / `G_electronic_p_v1` 保留为 rich-electronic sensitivity。

---

# 11. Electronic 金属/半导体二分不得人为修复

Electronic graph 中已有 metallic / semiconducting components。

禁止：
- 人工添加跨 component edge；
- 提高 k 直到强行连通；
- 加 epsilon bridge；
- 修改 Eg。

联合 graph 可通过 identity edge 间接连接，但 Electronic layer 内部结构不得改。

---

# 12. Phase P-B：Layer Strength Normalization

对每个 layer：

```text
strength_i = Σ_j W_ij
mean_strength_v = mean(strength_i)
W_scaled_v = W_v / mean_strength_v
```

目标：

```text
mean node strength ≈ 1
```

原 frozen graph 不覆盖，scaled graph 单独保存。

输出：

`data/audit/multiview_layer_scale.csv`

至少包含：

```text
view
N_nodes
N_edges
mean_strength_before
median_strength_before
mean_strength_after
scale_factor
```

第一版：

```text
alpha_structure = 1
alpha_Eg = 1
alpha_mass = 1
alpha_transport = 1
```

禁止用 PF 调权。

---

# 13. Phase P-C：建立 Partial Multilayer Graph

n-atlas 节点包括：

```text
Consensus nodes
Structure copies
Eg copies
electron-mass copies
n-Transport copies
```

p-atlas 对应替换为 hole mass 和 p-Transport。

跨 View 只允许：

```text
C_i ↔ V_i
```

权重：

`λ / m_i`

---

# 14. Phase P-D：λ 扫描

扫描：

```text
λ = 0.01
λ = 0.03
λ = 0.10
λ = 0.30
λ = 1.00
λ = 3.00
λ = 10.00
```

必要时在转变区间增加：

```text
0.05
0.20
0.50
2.00
```

禁止使用 PF 选择 λ。

---

# 15. Phase P-E：Joint Spectral Embedding

对每个 λ 构造 supra adjacency matrix `A`。

```text
D_ii = Σ_j A_ij
L_sym = I - D^(-1/2) A D^(-1/2)
```

求前 30–50 个 non-trivial eigenvectors。

至少保存：

```text
phi_1
phi_2
...
phi_20
```

二维只用于展示。

---

# 16. 每个 λ 必须同时评价 Alignment 与 Preservation

## Alignment

同一 JID 的 Consensus 和 View copy 是否靠近。

定义：

```text
T_i,v = distance(Consensus_i, View_i,v)
```

并归一化：

```text
T_norm(i,v)
=
T_i,v / median_joint_pair_distance
```

## Preservation

原冻结图中：

`NN_original_v(i)`

Joint coordinates 中、仅在同一 View copies 内：

`NN_joint_v(i)`

定义：

```text
P_i(v)
=
|NN_original_v(i) ∩ NN_joint_v(i)| / k
```

平均得到：

`P_v`

n-type 必须报告：

```text
P_structure
P_Eg
P_m_electron
P_n_transport
```

p-type 对应。

---

# 17. λ 选择规则

不根据二维图或 PF。

先排除任何使主要 View：

`P_v < 0.60`

的 λ。

优先寻找：

`P_v ≥ 0.70`

且 anchor tension 已进入下降平台的区域。

如果多个 λ 满足，选较小 λ。

必须报告完整 λ sensitivity。

若不存在满足条件的 λ：

`JOINT_MANIFOLD_SUPPORTED = False`

禁止继续调参直到得到成功。

---

# 18. Phase P-F：Structure Inclusion Test

最高优先级 ablation。

## n-Full

```text
Structure
+
Eg
+
electron mass
+
n-Transport
```

## n-Property

```text
Eg
+
electron mass
+
n-Transport
```

## p-Full

```text
Structure
+
Eg
+
hole mass
+
p-Transport
```

## p-Property

```text
Eg
+
hole mass
+
p-Transport
```

比较：
- Transport preservation；
- Eg preservation；
- Mass preservation；
- anchor tension；
- λ stable range；
- consensus neighborhood stability；
- spectral stability。

如果 Structure inclusion 明显破坏 property geometry，则 Structure 变为 linked auxiliary view。

不得为了课题原设想强行保留 Structure。

---

# 19. Phase P-G：决定 Unified 还是 Linked

必须明确输出：

```text
N_ATLAS_SUPPORTED = True / False
P_ATLAS_SUPPORTED = True / False
STRUCTURE_IN_COMMON_SPACE = True / False
```

如果 Supported，冻结：

```text
FINAL_N_ATLAS_V1
FINAL_P_ATLAS_V1
```

如果 Not Supported，冻结：

`FINAL_LINKED_MULTIVIEW_ATLAS_V1`

---

# 20. Phase Q：冻结 Joint Architecture

每个最终模型记录：

```text
model_name
carrier_type
included_views
alpha_per_view
lambda
layer_scale_factor
anchor_weight_rule
latent_dimension
graph_normalization
k
duplicate_policy
random_seed
```

保存 JSON。

---

# 21. Latent Dimension

不得预设为 2。

保存前 50 个 eigenvalues，检查：
- eigengap；
- spectral decay；
- participation ratio；
- 可选 TwoNN；
- 可选 MLE intrinsic dimension。

保存：

```text
d = 2
d = 3
d = 5
d = 10
d = 20
```

正式分析使用 validated dimension 或 d=10 baseline。

---

# 22. 最终 Atlas 使用 Consensus Node

最终：

`one JID = one Consensus coordinate`

同时保存各 View-copy coordinates，用于 tension。

对每个 JID 保存：

```text
view_count
available_views
transport_informed
mass_available
```

---

# 23. Core / Extended Atlas

## n-Core Atlas
仅包含具有真实 n-Transport 的约 806 个材料。

用于：
- PF mapping；
- transport landscape；
- tension；
- region analysis。

## n-Extended Atlas
显示全部 1103 JIDs。

没有 n-Transport 的材料必须：

`transport_informed = False`

p 型同样建立 Core / Extended。

---

# 24. Phase Q：Random Anchor Negative Control

保持所有单 View graph 不变，仅随机打乱 JID identity correspondence。

初步：

`200 permutations`

最终：

`1000 permutations`

比较：
- anchor tension；
- intra-view preservation；
- consensus-neighbor reproducibility；
- spectral structure；
- joint-space stability。

如果真实 JID 与 random 没有显著差异：

`JOINT_ALIGNMENT_VALID = False`

---

# 25. Phase Q：Layer Ablation

n-atlas 依次删除：

```text
Structure
Eg
electron mass
n-Transport
```

p-atlas 对应。

每次计算：

```text
consensus kNN overlap
distance rank correlation
Procrustes similarity
transport preservation
electronic preservation
anchor tension change
```

---

# 26. Duplicate Sensitivity

已有 `d_structure = 0` 的不同 JID。

主模型保留全部 JID。

建立：

`duplicate_group_id`

区分：

```text
exact duplicate
near duplicate
unique
```

再构建 exact-duplicate collapsed sensitivity dataset。

比较：
- consensus kNN overlap；
- distance rank correlation；
- PF landscape；
- tension rank；
- high-tension materials。

如果主要结论不变：

`DUPLICATE_SENSITIVITY_PASSED = True`

---

# 27. Phase R：正式 Joint Atlas

只有以下完成后才允许：
- λ 冻结；
- architecture 冻结；
- random-anchor passed；
- layer-ablation completed；
- duplicate sensitivity completed。

主图使用：

```text
Consensus Φ1
Consensus Φ2
```

n-type 依次着色：
- chemical family；
- Eg；
- electron mass；
- S_median；
- log_sigma_dom_geo；
- D_sigma；
- A_sigma_dom；
- JARVIS-defined PF；
- view_count；
- T_structure；
- T_transport。

p-type 对应。

---

# 28. PF External Validation

PF 必须在 architecture 与 λ 完全冻结后读取。

正式名称：

`JARVIS-defined PF`

禁止称：
- true directional PF；
- in-plane PF；
- x-PF；
- y-PF。

计算 Joint Atlas 中：

```text
Smoothness(log PF)
=
Σ W_ij × (logPF_i - logPF_j)^2
/
Σ W_ij
```

使用 1000 次 permutation null。

必须比较：

```text
Structure
Electronic
Transport
Joint
```

PF smoothness。

Joint 不要求一定优于 Transport。

---

# 29. 正式 Cross-View Tension

定义：

```text
T_i,v
=
distance(
Consensus_i,
ViewCopy_i,v
)
```

在最终高维 joint space 中计算。

输出：

```text
jid
formula
T_structure
T_Eg
T_mass
T_transport
T_mean
T_max
view_count
```

n：

`data/processed/joint_tension_n.csv`

p：

`data/processed/joint_tension_p.csv`

高 T_structure 只能解释为：

“Structure representation 与公共 electronic/transport representation 不一致。”

高 T_transport 可称：

“cross-view transport anomaly”。

---

# 30. Metal / Semiconductor 分析

定义：

```text
metal = Eg == 0
semiconductor = Eg > 0
```

回答：
1. Joint Atlas 是否仍形成 metal / semiconductor 区域？
2. Transport geometry 是否跨越边界？
3. high PF 是否主要集中在 semiconductor 区？
4. Seebeck sign-inconsistent materials 是否富集于 metallic / small-gap region？

---

# 31. Phase S：n/p Atlas 比较

两张 Atlas 分别训练，坐标轴不可直接相减。

共同 JID 上比较：

```text
consensus kNN overlap
pairwise distance Spearman correlation
local-neighborhood agreement
```

可使用 Procrustes 仅用于二维 figure 对齐。

---

# 32. Combined n+p Atlas

只有 n/p 主 Atlas 均完成后，才允许作为 sensitivity model。

Candidate Views：

```text
Structure
Eg
electron mass
hole mass
n-Transport
p-Transport
```

仍用 Consensus hub。

必须重新扫描 λ。

如果同时很好保留 n/p Transport：

`COMBINED_NP_ATLAS_SUPPORTED = True`

否则保持 n/p 两张独立 Atlas。

---

# 33. 本阶段禁止进入超晶格 Pair Ranking

即使已有：
- structure distance；
- electronic distance；
- transport distance；
- joint tension；

仍禁止最终 parent-pair ranking。

下一阶段还需引入：

```text
lattice compatibility
lattice mismatch
supercell commensurability
symmetry
interface orientation
band alignment
chemical compatibility
possible reconstruction
```

当前 Atlas 只提供 candidate parent pool。

---

# 34. Unified Atlas 成功标准

至少满足：

Criterion A：真实 JID anchors 显著优于 random anchors。

Criterion B：主要 View geometry 合理保留。

Criterion C：λ 存在稳定区间。

Criterion D：合理修改 k / representation 后 Consensus neighborhood 基本稳定。

Criterion E：duplicate collapse 不改变主要科学结论。

Criterion F：Structure inclusion / exclusion 的结论在合理参数范围内稳定。

---

# 35. Unified Atlas 失败条件

若出现：
- 没有合理 λ；
- random anchor 与真实 anchor 类似；
- Transport preservation 很低；
- Atlas 对 λ 极敏感；
- Structure inclusion 严重破坏 Electronic/Transport；
- duplicate collapse 导致整体重排；

则：

`UNIFIED_MANIFOLD_SUPPORTED = False`

并转为 Linked Multi-View Atlas。

---

# 36. 建议新增脚本

```text
scripts/
    31_build_mass_only_graphs.py
    32_scale_multiview_graphs.py
    33_build_consensus_hub_graph.py
    34_scan_anchor_lambda.py
    35_joint_feasibility_analysis.py
    36_structure_inclusion_test.py
    37_freeze_joint_architecture.py
    38_joint_spectral_embedding.py
    39_random_anchor_control.py
    40_layer_ablation.py
    41_duplicate_sensitivity.py
    42_joint_tension_analysis.py
    43_joint_external_pf_mapping.py
    44_compare_n_p_atlas.py
    45_combined_np_sensitivity.py
    46_phase_PS_summary.py
```

公共方法：

`scripts/multiview_utils.py`

---

# 37. 输出目录

```text
graphs/multiview/
├── n_atlas/
├── p_atlas/
└── combined_sensitivity/
```

Joint coordinates：

```text
manifolds/n_atlas_consensus.parquet
manifolds/n_atlas_all_nodes.parquet
manifolds/p_atlas_consensus.parquet
manifolds/p_atlas_all_nodes.parquet
```

---

# 38. λ Scan 输出

```text
data/audit/lambda_scan_n.csv
data/audit/lambda_scan_p.csv
```

必须包含：

```text
lambda
median_anchor_tension
P_structure
P_Eg
P_mass
P_transport
giant_component_fraction
spectral_gap
consensus_stability
```

---

# 39. Architecture Comparison

生成：

`data/audit/joint_architecture_comparison.csv`

至少包括：

```text
n-Full
n-Property
n-RichElectronic
p-Full
p-Property
p-RichElectronic
```

其他输出：

```text
data/audit/random_anchor_n.csv
data/audit/random_anchor_p.csv
data/audit/layer_ablation_n.csv
data/audit/layer_ablation_p.csv
data/audit/duplicate_sensitivity_n.csv
data/audit/duplicate_sensitivity_p.csv
data/audit/joint_pf_smoothness.csv
```

---

# 40. 必须单独生成报告

```text
reports/structure_inclusion_test.md
reports/joint_manifold_feasibility_n.md
reports/joint_manifold_feasibility_p.md
reports/multiview_layer_scaling.md
reports/lambda_selection.md
reports/random_anchor_control.md
reports/layer_ablation.md
reports/duplicate_sensitivity.md
reports/joint_tension_analysis.md
reports/joint_pf_analysis.md
reports/n_p_atlas_comparison.md
reports/phase_PS_summary.md
```

Joint feasibility 结论只能是：

`SUPPORTED`

或：

`NOT_SUPPORTED`

---

# 41. 必须生成的主要图

至少：

1. n-type λ vs anchor tension / preservation。
2. p-type λ vs anchor tension / preservation。
3. n-Full vs n-Property。
4. p-Full vs p-Property。
5. Real anchor vs Random anchor。
6. n-Atlas consensus map。
7. p-Atlas consensus map。
8. n-Atlas Eg coloring。
9. n-Atlas transport coloring。
10. n-Atlas PF coloring。
11. p-Atlas 对应图。
12. Structure tension map。
13. Transport tension map。
14. Layer ablation comparison。
15. n/p Atlas neighbor consistency。
16. Metal / semiconductor distribution。

图形规范：
- 白底；
- 不使用 3D；
- 坐标写 `Φ1`, `Φ2`；
- 不声称 Φ1/Φ2 是物理量；
- 缺失 property 灰色；
- transport-known 与 transport-missing 用不同 marker；
- colorbar 有单位；
- PF 写 `JARVIS-defined PF`；
- 同时保存 PNG 和 PDF/SVG。

---

# 42. 本阶段最终必须回答的问题

1. n-type joint manifold 是否 Supported？
2. p-type 是否 Supported？
3. n-type 最终 λ 是多少？
4. p-type 最终 λ 是多少？
5. λ 稳定区间是什么？
6. Structure 是否应该进入 n-atlas？
7. Structure 是否应该进入 p-atlas？
8. n-Full 和 n-Property 哪个更好？
9. p-Full 和 p-Property 哪个更好？
10. Transport preservation 最终是多少？
11. Eg preservation 是多少？
12. Mass preservation 是多少？
13. Structure preservation 是多少？
14. real anchors 相比 random anchors 显著性多大？
15. 哪个 layer ablation 影响最大？
16. duplicate collapse 是否影响结论？
17. n/p Atlas 是否真的属于不同公共 geometry？
18. Electronic–Transport 是否共享比 Structure–Transport 更强的共同几何？
19. metal / semiconductor boundary 是否仍存在？
20. PF 在 Joint Atlas 上是否显著平滑？
21. Joint PF smoothness 是否高于 Structure？
22. 是否高于 Electronic？
23. 是否高于 Transport？
24. 哪些材料 T_structure 最大？
25. 哪些材料 T_transport 最大？
26. high-tension materials 是否在参数变化下稳定？
27. Combined n+p Atlas 是否值得保留？
28. 最终应该称 Unified Atlas 还是 Linked Atlas？

---

# 43. 最终 Architecture Decision Table

生成：

| Model | Views | λ | Transport preservation | Electronic preservation | Structure preservation | Random-anchor | Decision |
|---|---|---:|---:|---:|---:|---|---|
| n-Full | Structure + Eg + m_e + Tn | ? | ? | ? | ? | ? | ? |
| n-Property | Eg + m_e + Tn | ? | ? | ? | — | ? | ? |
| n-Rich | Structure? + Electronic-n + Tn | ? | ? | ? | ? | ? | sensitivity |
| p-Full | Structure + Eg + m_h + Tp | ? | ? | ? | ? | ? | ? |
| p-Property | Eg + m_h + Tp | ? | ? | ? | — | ? | ? |
| p-Rich | Structure? + Electronic-p + Tp | ? | ? | ? | ? | ? | sensitivity |
| Combined | optional | ? | ? | ? | ? | ? | sensitivity |

所有问号必须来自真实运行结果。

---

# 44. 最终状态变量

必须写出：

```text
N_ATLAS_SUPPORTED = True / False
P_ATLAS_SUPPORTED = True / False
STRUCTURE_IN_N_COMMON_SPACE = True / False
STRUCTURE_IN_P_COMMON_SPACE = True / False
COMBINED_NP_ATLAS_SUPPORTED = True / False
UNIFIED_MANIFOLD_SUPPORTED = True / False
RANDOM_ANCHOR_CONTROL_PASSED = True / False
DUPLICATE_SENSITIVITY_PASSED = True / False
```

---

# 45. 最终 Architecture 输出格式

如果 n-Property 最优，例如：

```text
FINAL N-TYPE ATLAS V1

Consensus Identity Node
        |
        +--- Eg
        |
        +--- electron effective mass
        |
        +--- n-Transport

Structure:
linked auxiliary view

PF:
external label

kappa_e:
sensitivity only
```

这里只是示例，必须根据真实结果填写。

---

# 46. 如果 Unified Manifold 成立

冻结：

```text
FINAL_N_ATLAS_V1
FINAL_P_ATLAS_V1
```

下一阶段才做 Superlattice Parent Discovery。

如果不成立：

冻结：

`FINAL_LINKED_MULTIVIEW_ATLAS_V1`

下一阶段改为依据：

```text
Structure distance
Electronic distance
Transport distance
Cross-view disagreement
```

筛选材料对。

---

# 47. 本阶段禁止执行

禁止：
- ML 补数据；
- κ_L prediction；
- ZT 计算；
- 新 VASP；
- AMSET；
- BoltzTraP；
- phono3py；
- ShengBTE；
- superlattice generation；
- heterostructure generation；
- lattice matching；
- twist-angle search；
- interface optimization；
- 最终 parent-pair ranking。

---

# 48. STOP 条件

完成：

```text
mass-only graphs
layer scaling
consensus-hub graph
λ scan
joint feasibility
Structure inclusion test
random-anchor control
layer ablation
duplicate sensitivity
joint embedding
joint tension
PF external mapping
n/p comparison
optional combined sensitivity
```

以后：

**STOP**

不得进入超晶格设计。

---

# 49. 最终返回格式

## Executive Summary

### A. n-Atlas
`SUPPORTED / NOT SUPPORTED`

### B. p-Atlas
`SUPPORTED / NOT SUPPORTED`

### C. Structure
`IN COMMON SPACE / AUXILIARY ONLY`

### D. λ
给出：
```text
lambda_n
lambda_p
stable ranges
```

### E. Main Scientific Finding
回答：

> Electronic–Transport 是否形成比 Structure–Transport 更明显的共同材料几何？

## Final Architecture
给出最终结构图。

## Frozen Joint Models
列出：
```text
model
views
N
lambda
alpha
latent dimension
transport preservation
electronic preservation
structure preservation
random-anchor result
```

## PF External Validation
报告：
```text
PF smoothness on Structure
PF smoothness on Electronic
PF smoothness on Transport
PF smoothness on Joint
```

## Cross-View Physics
重点解释：
```text
Structure ↔ Electronic
Structure ↔ Transport
Electronic ↔ Transport
n ↔ p
```

## Stable Anomalous Materials
列出：
- high Structure tension；
- high Transport tension；
- 参数稳定异常点。

只能称 anomalous materials，不得称 superlattice candidates。

## Final Status Variables
完整列出所有状态变量。

## Ready for Superlattice Parent Discovery?
最后只回答：
`YES`
或
`NO`

如果 YES，下一阶段才进入：

**structure-compatible + electronic-compatible + transport-contrasting parent pair discovery**

---

# 50. 本阶段最终科学原则

当前最重要的不是：

> 如何把所有材料强行塞进一张漂亮二维图。

而是：

> 是否存在一个可以同时尊重不同物理 View 原始局部几何的公共材料表示。

Structure 与 Transport 已经被实证为几乎独立。

因此 Structure 是否进入最终公共空间，必须由 preservation / alignment 结果决定。

如果 Structure 最终被排除，这不是失败。

如果统一流形最终不成立，也不是失败。

有可能真正的科学结论正是：

> 二维材料的结构相似性、电子相似性和热电输运相似性对应不同但彼此关联的材料空间，而不是一个单一全局流形。

因此：

**真实数据优先于预设模型。**

现在开始执行：

```text
Phase P
Partial Multi-View Alignment
        ↓
Phase Q
Joint Feasibility Validation
        ↓
Phase R
Unified / Linked Atlas
        ↓
Phase S
PF External Mapping + n/p Comparison
        ↓
STOP
```

完成 Phase S 后停止，并按规定格式返回全部结果。
