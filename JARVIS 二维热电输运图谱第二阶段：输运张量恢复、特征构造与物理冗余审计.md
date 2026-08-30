# JARVIS 二维热电输运图谱第二阶段：输运张量恢复、特征构造与物理冗余审计

## 一、当前研究状态

第一阶段已经完成。

当前数据库来自：

JARVIS `dft_2d`

通过官方 JARVIS OPTIMADE API 获取。

当前材料总数：

`N_total = 1103`

已经确认的主要 property coverage：

| Property | Available | Coverage |
|---|---:|---:|
| Structure | 1103 | 100% |
| OptB88vdW band gap | 1103 | 100% |
| Effective mass | 678 | 61.5% |
| n-Seebeck | 806 | 73.1% |
| p-Seebeck | 802 | 72.7% |
| n-PF | 802 | 72.7% |
| p-PF | 802 | 72.7% |
| n-conductivity | 802 | 72.7% |
| p-conductivity | 802 | 72.7% |
| n-kappa_e | 802 | 72.7% |
| p-kappa_e | 802 | 72.7% |
| kappa_L | 0 | 0% |
| ZT | 0 | 0% |

因此当前研究正式定义为：

**2D Thermoelectric Transport Atlas**

而不是 Full-ZT Atlas。

当前阶段不研究 kappa_L 和 ZT。

---

# 二、已经确认的 JARVIS BoltzTraP 数据定义

通过检查：

`jarvis/db/vasp_to_xml.py`

中的：

`boltztrap_data()`

已经确认以下事实。

所有热电输运数据均来自固定条件：

`T = 600 K`

`|carrier concentration| = 1e20 cm^-3`

其中：

- n 型对应电子掺杂；
- p 型对应空穴掺杂。

不是温度平均。

不是 carrier-concentration 平均。

不是温度扫描 response surface。

---

## 2.1 Seebeck

原始 BoltzTraP tensor 为：

`3 x 3`

JARVIS 对其执行：

`eigvals(seeb_tensor * 1e6)`

得到 3 个主值。

单位：

`microV/K`

即：

`µV/K`

---

## 2.2 Conductivity

原始 conductivity tensor 为：

`3 x 3`

数据库处理：

`eigvals(cond_tensor) / 1e14`

得到 3 个主值。

与固定 relaxation-time convention 有关。

当前不要把该值解释成严格实验 conductivity。

但由于所有材料使用统一协议：

**允许进行数据库内部相对比较、排序和相似性分析。**

---

## 2.3 Power Factor

源码计算逻辑：

`PF_i = S_i^2 * sigma_i / 1e6`

即针对对应的 principal values 计算。

PF 因此不是独立原始物理自由度，而是由：

`S`

和：

`sigma`

派生。

因此主流形中默认：

`PF = external validation property`

不要让 PF 与 S 和 sigma 同时作为等权输入。

---

## 2.4 Electronic thermal conductivity

原始：

`kappa_e`

同样经过 3 x 3 tensor diagonalization。

但 JARVIS 对其没有执行与 conductivity 相同的 `/1e14` 缩放。

因此当前：

`kappa_e`

的绝对数值不能直接解释成标准 W/(m K)。

它仍然可以在统一计算协议内部：

- 比较材料；
- 计算排名；
- 建立相似性；
- 分析相关性。

---

# 三、当前发现的一个重要数据损失问题

JARVIS 原始 FigShare JSON 中：

`S`

`sigma`

`PF`

`kappa_e`

理论上分别保存 3 个 tensor eigenvalues。

但是当前通过 OPTIMADE API 获取的数据中，这 3 个值已经被平均为一个 scalar。

例如：

`S_mean = (S1 + S2 + S3) / 3`

因此当前 OPTIMADE 数据已经丢失：

- principal-value distribution；
- anisotropy；
- max/min；
- tensor spectral spread。

这对于二维材料尤其重要。

因此本轮第一优先级是：

**尝试从同一个 JARVIS 数据体系恢复原始 3 eigenvalues。**

注意：

这不属于使用其他数据库补数据。

它仍然属于恢复同一 JARVIS 记录的完整 representation。

---

# 四、本轮严格目标

本轮只完成：

### Phase G

尝试恢复 JARVIS 原始 3 个输运本征值。

### Phase H

完成输运字段内部结构和 JID 交集审计。

### Phase I

建立物理一致、排列不变的 transport descriptors。

### Phase J

完成相关性、冗余性和 anisotropy 审计。

### Phase K

确定下一阶段真正需要建立哪些 physical views。

完成 Phase K 后：

**STOP**

本轮不要开始：

- SOAP manifold；
- UMAP；
- Diffusion Map；
- supra graph；
- joint spectral embedding；
- candidate pair selection。

先把输运 representation 搞清楚。

---

# 五、Phase G：尝试恢复原始 3 eigenvalues

## G1. 不要继续使用 OPTIMADE scalar 作为唯一来源

当前 OPTIMADE scalar 数据保留。

它作为：

`mean-value fallback dataset`

不要删除。

但是首先尝试寻找同一个 JID 对应的原始三本征值。

---

# 六、Phase G1：检查 JARVIS 官方 Python 接口和缓存机制

检查当前安装的：

`jarvis-tools`

源码。

重点检查：

`jarvis/db/figshare.py`

以及：

`get_db_info()`

确认：

`dft_2d`

的：

- FigShare file ID；
- archive URL；
- JSON filename；
- local cache filename。

记录到：

`reports/raw_transport_source_audit.md`

---

# 七、Phase G2：检查当前机器是否已经存在 jarvis-tools 缓存

搜索：

- 用户 home；
- conda environment；
- jarvis cache；
- 当前项目；
- Python site-packages；
- `/tmp`；
- `$HOME/.jarvis`；
- `$HOME/.cache`；

是否已经存在：

`d2-12-12-2022.json`

或者对应 ZIP。

如果存在：

不要重新下载。

首先：

- 计算 SHA256；
- 检查文件完整性；
- 检查 JSON；
- 确认材料数量；
- 确认是否真的包含三个本征值。

---

# 八、Phase G3：尝试其他 JARVIS 官方接口

如果 FigShare 当前环境仍然 403：

只允许尝试 JARVIS 官方来源。

例如：

- JARVIS REST API；
- JARVIS-DFT backend；
- JARVIS OPTIMADE additional fields；
- 官方 static files；
- 官方 GitHub 中公开的数据入口；
- NIST-hosted mirrors。

禁止使用：

- Kaggle 镜像；
- 不明 GitHub fork；
- 网盘；
- 第三方转载数据库。

对每一个候选接口记录：

| Source | Official? | Accessible? | Returns eigenvalues? |
|---|---|---|---|

输出：

`reports/jarvis_transport_source_probe.csv`

---

# 九、Phase G4：检查 OPTIMADE 是否隐藏保留数组字段

不要只查看当前导出的 scalar property。

检查每个 structure entry 的：

- attributes；
- provider-specific fields；
- `_jarvis_*` fields；
- raw JSON response；

搜索关键词：

`seeb`

`cond`

`power`

`kappa`

`boltz`

`tensor`

`eigen`

`transport`

检查是否存在：

- array field；
- nested dictionary；
- raw-value field；
- tensor field；

只是当前转换脚本没有读取。

如果存在：

优先使用原始数组。

---

# 十、Phase G5：如果成功获得原始 3 eigenvalues

必须针对每个 JID 保存：

`n_S_eigs`

`p_S_eigs`

`n_sigma_eigs`

`p_sigma_eigs`

`n_PF_eigs`

`p_PF_eigs`

`n_kappa_e_eigs`

`p_kappa_e_eigs`

每个字段应严格包含：

`3 values`

---

# 十一、不要假设 eigvals 的顺序有固定方向意义

这是非常重要的物理和数学问题。

JARVIS 使用：

`np.linalg.eigvals()`

因此原始三个值只是 tensor eigenvalues。

它们不是：

`xx`

`yy`

`zz`

也不是固定：

`x`

`y`

`z`

方向。

不同材料之间：

`eigvalue #1`

不保证对应同一个 physical axis。

因此禁止直接使用：

`[eig1, eig2, eig3]`

作为跨材料 feature vector，

除非首先执行明确、统一的排序规则。

---

# 十二、首先检查 eigenvalues 是否始终为实数

对于所有 tensor eigenvalues：

检查：

`imaginary part`

统计：

- maximum absolute imaginary component；
- number of records with non-negligible imaginary component。

建议容差：

`abs(Im) < 1e-8`

如果均满足：

转换为 real。

如果存在明显 complex eigenvalues：

停止。

必须检查 tensor 是否非对称，或者数据处理是否有问题。

---

# 十三、建立 permutation-invariant tensor descriptors

推荐不要把三个本征值简单按数组顺序使用。

对 eigenvalue set：

`{v1, v2, v3}`

构造以下 permutation-invariant descriptors：

`mean`

`median`

`std`

`min`

`max`

`range = max - min`

`abs_mean`

`abs_max`

`RMS`

此外根据 property 的物理含义构造 anisotropy descriptor。

---

# 十四、Conductivity 与 kappa_e 的 anisotropy

因为：

`sigma > 0`

理论上 conductivity principal values 应为正。

定义：

`anisotropy_ratio = max / min`

但为了防止极小值：

只有当：

`min > epsilon`

时使用。

否则设置：

`anisotropy_ratio = NaN`

并记录原因。

也可以使用更稳定指标：

`anisotropy_log = log10(max / min)`

---

# 十五、Seebeck 的 anisotropy 不能简单使用 max/min

因为 Seebeck 可以为负。

对 Seebeck 推荐：

`S_mean`

`S_std`

`S_min`

`S_max`

`S_range`

`abs_S_mean`

`abs_S_max`

以及：

`S_relative_spread = std / (abs(mean) + epsilon)`

不要直接使用：

`max / min`

作为主要 Seebeck anisotropy。

---

# 十六、PF 的处理

PF 是派生性质。

保存其：

- mean；
- std；
- min；
- max；
- anisotropy；

用于描述和 validation。

但是默认不要让 PF 进入主 transport graph。

---

# 十七、如果无法恢复 3 eigenvalues

如果所有 JARVIS 官方接口均无法获得原始数组：

不要中止整个项目。

正式记录：

`RAW_EIGENVALUES_AVAILABLE = False`

继续使用当前 OPTIMADE scalar mean。

但是所有后续文档必须明确写：

**Transport properties represent the arithmetic mean of three principal tensor eigenvalues.**

此时：

禁止研究：

- transport anisotropy；
- tensor principal-value dispersion；
- orientation dependence。

对应章节全部：

`SKIPPED_RAW_TENSOR_UNAVAILABLE`

---

# 十八、验证 OPTIMADE scalar 是否真的是 eigenvalue mean

如果成功恢复原始三本征值：

必须随机抽取至少：

`50 materials`

验证：

`OPTIMADE_scalar`

是否等于：

`mean(raw_3_eigenvalues)`

记录误差：

`absolute_error`

`relative_error`

输出：

`reports/optimade_mean_validation.csv`

并报告：

- MAE；
- max error；
- number of exact matches；
- number of mismatches。

只有通过以后，才能确认当前 scalar 的数学定义。

---

# 十九、Phase H：JID 交集审计

分别建立以下 JID sets：

`D_S_n`

`D_S_p`

`D_sigma_n`

`D_sigma_p`

`D_kappa_e_n`

`D_kappa_e_p`

`D_PF_n`

`D_PF_p`

`D_mstar`

`D_Eg`

---

# 二十、计算所有主要 view 的交集数量

至少计算：

`|D_S_n ∩ D_sigma_n|`

`|D_S_n ∩ D_kappa_e_n|`

`|D_sigma_n ∩ D_kappa_e_n|`

`|D_S_n ∩ D_sigma_n ∩ D_kappa_e_n|`

以及 p 型对应结果。

还要计算：

`|D_mstar ∩ D_transport_n|`

`|D_mstar ∩ D_transport_p|`

`|D_Eg ∩ D_transport_n|`

`|D_Eg ∩ D_transport_p|`

其中：

`D_transport_n = D_S_n ∩ D_sigma_n ∩ D_kappa_e_n`

`D_transport_p = D_S_p ∩ D_sigma_p ∩ D_kappa_e_p`

---

# 二十一、建立 JID overlap matrix

输出：

`data/audit/view_jid_overlap.csv`

形式：

| View A | View B | N_A | N_B | N_intersection | Fraction_A | Fraction_B |
|---|---|---:|---:|---:|---:|---:|

同时生成 heatmap。

这一步用于确认：

虽然多个字段都有约 802 条数据，

它们是否真的是同一批 JIDs。

不要根据数量相同推断集合相同。

---

# 二十二、Phase H2：字段数值质量检查

对每一个 transport property：

分别对 n 和 p：

统计：

- count；
- NaN；
- Inf；
- zero；
- positive；
- negative；
- min；
- 1st percentile；
- median；
- 99th percentile；
- max。

输出：

`data/audit/transport_numeric_audit.csv`

---

# 二十三、Seebeck sign 检查

检查：

n-type Seebeck

是否大多数为负。

检查：

p-type Seebeck

是否大多数为正。

统计 sign violation 数量。

不要自动删除 sign violation。

输出材料 JID，

后续人工检查。

这些异常可能来源于：

- band complexity；
- convention；
- numerical issue；
- carrier definition。

---

# 二十四、Conductivity 正值检查

理论上 conductivity principal values 应：

`>= 0`

检查：

- negative values；
- zero values；
- extremely small values。

出现负值必须单独报告。

---

# 二十五、kappa_e 正值检查

同样检查：

`kappa_e >= 0`

异常值不要静默删除。

---

# 二十六、PF 正值检查

因为 PF 来源于：

`S^2 * sigma`

理论上：

`PF >= 0`

出现负值表示数据处理存在异常。

---

# 二十七、Phase I：建立输运特征表

如果成功恢复 3 eigenvalues：

分别建立：

`features/transport/n_transport_tensor_features.parquet`

`features/transport/p_transport_tensor_features.parquet`

每个 JID 包括：

### Seebeck

`S_mean`

`S_std`

`S_min`

`S_max`

`S_range`

`S_abs_mean`

`S_relative_spread`

### Conductivity

`sigma_mean`

`sigma_std`

`sigma_min`

`sigma_max`

`sigma_anisotropy_log`

### Electronic thermal conductivity

`kappa_e_mean`

`kappa_e_std`

`kappa_e_min`

`kappa_e_max`

`kappa_e_anisotropy_log`

### PF

仅保存用于 validation：

`PF_mean`

`PF_std`

`PF_min`

`PF_max`

`PF_anisotropy_log`

---

# 二十八、如果只有 OPTIMADE scalar

则分别建立：

`n_transport_scalar_features.parquet`

和：

`p_transport_scalar_features.parquet`

只包括：

`S_mean`

`sigma_mean`

`kappa_e_mean`

以及 external label：

`PF_mean`

此时不要制造不存在的 anisotropy features。

---

# 二十九、不要直接使用原始 sigma 和 kappa_e 进行欧氏距离

因为它们可能跨越多个数量级。

对：

`sigma`

使用：

`log10(sigma + epsilon_sigma)`

对：

`kappa_e`

使用：

`log10(kappa_e + epsilon_kappa)`

epsilon 必须根据数据尺度定义并记录。

不要随意使用：

`1e-10`

而不解释。

优先根据最小正常正值确定，例如：

`epsilon = min_positive / 10`

并做 sensitivity check。

---

# 三十、Seebeck 默认不要取绝对值

Seebeck 的正负包含 carrier-type 信息。

由于 n 型和 p 型已经分别分析：

n-type 保留原始 S_n。

p-type 保留原始 S_p。

只在 anisotropy 描述时使用 abs-based statistics。

---

# 三十一、Phase J：相关性与冗余审计

这是开始真正流形分析之前最关键的一步之一。

分别对 n-type 和 p-type 计算：

Pearson correlation。

Spearman correlation。

重点分析：

`S`

`sigma`

`kappa_e`

`PF`

---

# 三十二、首先检查 sigma 和 kappa_e

特别关注：

`corr(log_sigma, log_kappa_e)`

如果：

`|Spearman rho| > 0.95`

或者极度接近 1，

说明二者在当前固定 T 和 doping 下提供的信息高度重复。

此时不要让：

`sigma`

和：

`kappa_e`

作为两个完全独立、等权的 property layers。

---

# 三十三、检查 PF 的派生冗余

计算：

`corr(PF, S)`

`corr(PF, log_sigma)`

并验证数据库 PF 与：

`S^2 * sigma`

的数值关系。

注意：

如果使用 mean eigenvalues：

一般不能假设：

`mean(PF_i) = mean(S_i)^2 * mean(sigma_i)`

因为：

`mean(S_i^2 * sigma_i)`

通常不等于：

`mean(S_i)^2 * mean(sigma_i)`

所以不要用 OPTIMADE mean S 和 mean sigma 重建 mean PF。

如果恢复原始对应 eigenvalues，

才可以逐 principal value 验证源码公式。

---

# 三十四、PF 默认作为 external performance label

主模型：

`Transport View A`

仅使用：

`S`

`sigma`

`kappa_e`

或其去冗余版本。

PF 不参与 embedding。

然后在建立 transport space 以后：

按 PF 着色。

核心问题：

**高 PF 是否自然聚集在由基础输运量构建的空间中？**

这样避免把目标性质直接输入后又“发现”高 PF 区域。

---

# 三十五、建立 PF sensitivity model

后续可以建立：

`Transport View B`

包含：

`S`

`sigma`

`kappa_e`

`PF`

仅作为 sensitivity analysis。

比较 View A 与 View B：

- neighbor overlap；
- graph distance rank correlation；
- embedding stability。

如果加入 PF 导致空间大幅变化，

说明 PF 在重复或主导已有输运信息。

主论文优先使用 View A。

---

# 三十六、Phase J2：Principal Component / Effective Rank 审计

在不进行最终 manifold 建模的情况下，

允许对标准化 transport features 做：

PCA diagnostic。

目的仅仅是检查 feature redundancy。

报告：

- explained variance ratio；
- cumulative explained variance；
- effective rank。

如果：

第一主成分已经解释 >90% 信息，

说明当前 transport representation 接近一维。

此时不要夸大称为复杂 transport manifold。

如果需要多个主成分，

说明确实存在多维 transport geometry。

PCA 这里只是 audit 工具，

不是最终材料图。

---

# 三十七、Phase J3：各向异性信息价值

仅当恢复原始 3 eigenvalues 时执行。

比较：

### Mean-only representation

例如：

`[S_mean, log_sigma_mean, log_kappa_mean]`

与：

### Tensor-spectrum representation

例如：

`[S_mean, S_std, S_range, log_sigma_mean, sigma_anisotropy_log, log_kappa_mean, kappa_anisotropy_log]`

计算：

- PCA effective dimension；
- pairwise-distance rank correlation；
- kNN overlap。

回答：

**保留 tensor eigenvalue spectrum 是否真正改变材料邻域？**

如果改变很小，

mean 已足够。

如果改变明显，

说明 anisotropy 是重要输运自由度。

---

# 三十八、不要把 3 个 eigenvalues 当成空间方向

即使成功获得：

`v1, v2, v3`

也禁止解释：

`v1 = x`

`v2 = y`

`v3 = z`

因为当前数据库只保存 eigenvalues，

没有保存对应 eigenvectors。

因此只能研究：

**principal-value spectrum**

不能研究：

**crystallographic transport direction**

---

# 三十九、二维材料的一个特殊问题

对于真正二维材料，

理论上通常存在：

- 两个 in-plane principal directions；
- 一个 out-of-plane-like direction。

但是仅凭 3 eigenvalues，

无法可靠知道哪个 eigenvalue 对应 out-of-plane，

除非同时获得 eigenvectors 或原始 tensor 与晶格坐标系。

因此：

**禁止简单删除三个 eigenvalues 中最小的一个并称其为 z 方向。**

只有获得 eigenvectors 或原始 3 x 3 tensor 后，

才能进行真正的 in-plane/out-of-plane 分离。

---

# 四十、如果可以获得原始 3 x 3 tensor

如果 JARVIS 官方接口能够提供原始 tensor：

优先级高于 eigenvalues。

此时保存：

`3 x 3 Seebeck tensor`

`3 x 3 conductivity tensor`

`3 x 3 kappa_e tensor`

同时保存晶格矩阵。

但是本轮先只进行数据审计。

不要立刻展开完整方向性输运分析。

在报告中标记：

`RAW_TENSOR_AVAILABLE = True`

后续再单独设计 2D tensor projection 方法。

---

# 四十一、Phase K：确定最终 Physical Views

完成以上审计后，

不要根据原计划机械建立很多 layers。

根据真实数据和冗余结果选择。

第一候选架构：

### View 1：Structure

覆盖：

`1103`

未来使用：

SOAP + composition + local geometry。

### View 2：Electronic

主要使用：

`OptB88vdW gap`

加：

`effective mass`

仅在有效质量存在的子集中形成 richer electronic view。

### View 3：n-type Transport

候选：

`S_n`

`log_sigma_n`

`log_kappa_e_n`

以及 tensor spectral descriptors，如果能够恢复。

不输入 PF_n。

### View 4：p-type Transport

候选：

`S_p`

`log_sigma_p`

`log_kappa_e_p`

以及 tensor spectral descriptors。

不输入 PF_p。

---

# 四十二、如果 sigma 与 kappa_e 高度冗余

如果审计发现：

`|rho(log_sigma, log_kappa_e)| >= 0.95`

则建立两个候选模型。

### Model T1

`[S, log_sigma]`

### Model T2

`[S, log_kappa_e]`

比较：

- pairwise distance；
- kNN overlap；
- PF smoothness。

选取更简洁且更稳定的 representation。

不要重复计权同一个输运自由度。

---

# 四十三、Electronic View 也要注意重复信息

如果 effective mass 本身包含：

- electron mass；
- hole mass；
- multiple directions；

不要简单把全部 flatten 后使用。

先检查结构。

同样建立：

- permutation-invariant descriptors；
- n/p 对应关系；
- dimensionality audit。

如果 effective mass 数据与当前 n/p transport carrier type 可以合理对应，

再进入 Electronic View。

---

# 四十四、本轮必须输出的数据文件

至少生成：

`data/audit/view_jid_overlap.csv`

`data/audit/transport_numeric_audit.csv`

`data/audit/transport_correlation_pearson.csv`

`data/audit/transport_correlation_spearman.csv`

`data/audit/transport_pca_diagnostics.csv`

`reports/jarvis_transport_source_probe.csv`

`reports/raw_transport_source_audit.md`

`reports/transport_redundancy_analysis.md`

`reports/round2_summary.md`

---

# 四十五、如果恢复原始本征值，还必须生成

`data/processed/transport_eigenvalues.parquet`

以及：

`features/transport/n_transport_tensor_features.parquet`

`features/transport/p_transport_tensor_features.parquet`

另外生成：

`reports/optimade_mean_validation.csv`

---

# 四十六、如果无法恢复原始本征值

则明确生成：

`reports/raw_tensor_recovery_failed.md`

内容说明：

- 尝试了哪些官方接口；
- HTTP 状态；
- 为什么失败；
- 为什么最终继续使用 OPTIMADE mean；
- 后续分析有哪些限制。

不要将这视为程序失败。

这只是决定 representation 层级。

---

# 四十七、本轮最终必须回答的问题

完成以后，必须明确回答以下问题。

1. 是否成功恢复原始 3 eigenvalues？

2. 是否能够获取原始 3 x 3 tensor？

3. OPTIMADE scalar 是否严格等于 3 eigenvalues 的 arithmetic mean？

4. n-S、n-sigma、n-kappa_e 是否对应同一批 JID？

5. p 型对应集合是否相同？

6. sigma 与 kappa_e 的相关系数是多少？

7. PF 与基础输运量的冗余程度多大？

8. Transport View 的 effective dimension 大约是多少？

9. tensor spectral anisotropy 是否明显改变材料邻域？

10. 最终推荐使用哪些 transport variables？

11. 最终推荐建立几个 physical views？

12. PF 是否继续作为 external performance label？

---

# 四十八、本轮最终必须生成决策表

生成类似：

| Candidate View | Features | N | Keep? | Reason |
|---|---|---:|---|---|
| Structure | structure descriptors | 1103 | YES | backbone |
| Electronic | Eg + effective mass | ? | ? | ? |
| n-Transport | S_n + log sigma_n | ? | ? | ? |
| n-Transport extended | S_n + log sigma_n + log kappa_e_n | ? | ? | redundancy-dependent |
| p-Transport | S_p + log sigma_p | ? | ? | ? |
| PF-n | PF_n | ? | NO as input | external validation |
| PF-p | PF_p | ? | NO as input | external validation |
| kappa_L | unavailable | 0 | NO | not available |
| ZT | unavailable | 0 | NO | not available |

所有问号必须来自真实审计。

---

# 四十九、停止条件

本轮完成：

Phase G

Phase H

Phase I

Phase J

Phase K

后：

**STOP**

不要开始：

- SOAP 计算；
- kNN graph；
- UMAP；
- Diffusion Map；
- partial manifold alignment；
- supra adjacency；
- joint Laplacian；
- unified atlas。

必须先把：

**最终 physical views 和 transport representation**

确认下来。

---

# 五十、下一阶段预告，但当前不要执行

等这一轮结果确认以后，

下一阶段才进入：

### Phase L

Structure descriptors。

### Phase M

Structure similarity graph。

### Phase N

Electronic、n-Transport、p-Transport similarity graphs。

### Phase O

Single-view geometry validation。

### Phase P

JID-anchor multilayer alignment。

### Phase Q

Joint spectral / diffusion embedding。

### Phase R

Unified 2D Thermoelectric Transport Atlas。

### Phase S

PF external mapping。

### Phase T

Structure-close / transport-far candidate pairs。

---

# 五十一、当前最重要的方法原则

不要追求 property layer 越多越好。

应该追求：

**彼此独立且物理意义清楚的 view。**

如果：

sigma 和 kappa_e 几乎完全共线，

它们不应该因为名字不同就被算作两个独立的物理空间。

如果：

PF 已经由 S 和 sigma 派生，

它不应该在主输运流形中再次获得与基础量相同的权重。

如果：

原始 tensor eigenvalues 可以恢复，

优先保留其各向异性信息。

如果只能获取均值，

诚实地把当前 representation 定义为：

**mean transport property space**

而不是 tensor-resolved manifold。

最终希望建立的是：

**Structure View + Electronic View + n-Transport View + p-Transport View**

通过 JID identity anchors 映射到统一材料空间。

PF 作为外部性能标签验证：

> 基础输运空间中是否自然形成 high-PF region。

当前立即从：

**Phase G：恢复和审计原始 JARVIS 输运 tensor information**

开始。

完成 Phase K 后 STOP，并返回完整审计报告。