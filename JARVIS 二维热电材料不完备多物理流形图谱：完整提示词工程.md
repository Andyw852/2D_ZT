# JARVIS 二维热电材料不完备多物理流形图谱：完整提示词工程

## 1. 研究任务

你现在是我的材料信息学、二维热电材料、数据工程、图机器学习和流形学习研究助手。

我要开展一个完整、可复现、具有论文发表潜力的研究项目：

**利用 JARVIS 数据库中已有的二维材料数据，在不进行新的第一性原理计算、不预测缺失性质、不插值缺失数据、不利用其他数据库填补缺失值的前提下，分别建立结构、电子结构和热电输运性质对应的材料相似性空间，并通过不完备多视图流形对齐或多层图嵌入，将这些覆盖范围不同的物理视图映射到一个统一的二维材料公共潜在空间。**

最终目标是建立：

**2D Thermoelectric Materials Atlas**

中文名称：

**二维热电材料多物理统一图谱**

当前阶段的核心不是预测新的材料性质，也不是强行获得完整 ZT，而是研究：

> 同一种二维材料在结构空间、电子空间和热电输运空间中的邻域关系是否一致，以及这些不同物理空间之间能否通过材料身份进行统一对齐。

---

# 2. 核心科学问题

需要系统回答以下问题。

1. 二维材料是否存在稳定的结构和成分低维流形？

2. 结构相似的二维材料是否具有相似的电子结构？

3. 结构相似是否意味着 Seebeck 系数相似？

4. Seebeck 相似是否意味着功率因子 PF 相似？

5. 不同热电性质是否对应不同的材料局部邻域？

6. 不同性质具有不同的数据覆盖率时，能否利用共同的材料 ID 作为锚点进行流形对齐？

7. 能否在完全不补全缺失性质的情况下建立统一材料潜在空间？

8. 高 Seebeck、高 PF、合适带隙和合适有效质量区域，在公共材料空间中是重合、相邻还是彼此竞争？

9. 是否存在：

`d_struct(A,B) ≪ 1`

但：

`d_transport(A,B) ≫ 1`

的材料对？

也就是说：

**两个材料结构非常相似，但热电输运行为显著不同。**

10. 这些“结构相似、输运不同”的材料对，是否能够成为后续二维超晶格或异质结构设计的重要候选？

---

# 3. 第一阶段只使用一个数据库

当前阶段只允许使用：

`JARVIS curated dft_2d`

暂时禁止混入：

- C2DB
- Materials Project
- 2DMatPedia
- Alexandria
- 自己的 VASP 数据
- 文献人工整理数据
- 实验数据
- ML 预测数据

这样做的原因是：

**尽量保持同一数据库和统一计算协议，避免数据库之间的方法差异被误认为材料本身的物理差异。**

---

# 4. 不允许补充或预测缺失数据

如果一个材料存在：

- structure
- band gap
- Seebeck

但没有：

- κ_L

则必须保留：

`κ_L = missing`

禁止使用：

- mean imputation
- median imputation
- KNN imputation
- zero filling
- matrix completion
- regression prediction
- machine-learning prediction
- surrogate model

核心原则：

**有该性质，就在对应性质层中建立节点。**

**没有该性质，就不存在这个性质节点。**

缺失值绝不能人为设置成 0。

---

# 5. 为什么不能只使用所有性质都完整的材料

假设数据库总共有 N 个二维材料。

结构数据集合记为：

`D_struct = {全部二维材料}`

具有 Seebeck 数据的材料集合：

`D_S ⊆ D_struct`

具有 PF 数据的材料集合：

`D_PF ⊆ D_struct`

具有有效质量数据的材料集合：

`D_m* ⊆ D_struct`

其他物性也是如此。

禁止简单取：

`D_complete = D_struct ∩ D_S ∩ D_PF ∩ D_m* ∩ ...`

然后只研究这些完整材料。

原因是这样会：

- 大量减少数据；
- 丢失数据库已有信息；
- 产生 complete-case selection bias；
- 使 κ_L 等低覆盖率性质成为整个数据规模的瓶颈。

正确方法是：

**Incomplete Multi-View Learning**

或：

**Partial Multi-View Manifold Alignment**

即：

不同物理视图允许拥有不同数量的材料。

---

# 6. 应该下载哪个 JARVIS 数据集

主要使用：

```python
data("dft_2d")
```

不要把：

```python
data("dft_2d_2021")
```

作为主数据库。

`dft_2d_2021` 后续可以用于版本稳定性验证。

当前也不要下载：

```python
data("dft_3d")
```

因为当前研究对象仅为二维材料。

---

# 7. 建立项目目录

创建：

```text
jarvis_2d_te_atlas/
│
├── README.md
├── environment.yml
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   └── jarvis/
│   ├── audit/
│   ├── interim/
│   └── processed/
│
├── structures/
│
├── features/
│   ├── structure/
│   ├── electronic/
│   └── transport/
│
├── graphs/
├── manifolds/
├── figures/
├── reports/
│
└── scripts/
    ├── 00_environment.py
    ├── 01_verify_database.py
    ├── 02_download_dft2d.py
    ├── 03_validate_archive.py
    ├── 04_schema_audit.py
    ├── 05_property_coverage.py
    ├── 06_structure_validation.py
    ├── 07_build_structure_features.py
    ├── 08_build_property_views.py
    ├── 09_build_single_view_graphs.py
    ├── 10_single_view_embeddings.py
    ├── 11_build_supra_graph.py
    ├── 12_joint_spectral_embedding.py
    ├── 13_validate_alignment.py
    ├── 14_property_mapping.py
    ├── 15_neighbor_analysis.py
    └── 16_final_report.py
```

---

# 8. Step 0：建立 Python 环境

推荐 Python：

`Python 3.11`

创建 Conda 环境：

```bash
conda create -n te_manifold python=3.11 -y
conda activate te_manifold
```

安装：

```bash
python -m pip install -U jarvis-tools numpy pandas scipy scikit-learn matplotlib pyarrow pymatgen ase networkx tqdm dscribe umap-learn requests
```

然后：

```bash
mkdir -p reports
python --version > reports/python_version.txt
pip freeze > reports/pip_freeze.txt
```

所有后续程序必须能够在这个环境中重新执行。

---

# 9. Step 1：验证 JARVIS 当前真实数据源

禁止把 FigShare 下载地址直接硬编码到程序中。

首先读取当前安装版本 `jarvis-tools` 的官方配置：

```python
from jarvis.db.figshare import get_db_info

info = get_db_info()

if "dft_2d" not in info:
    raise RuntimeError(
        "Current jarvis-tools does not provide dft_2d."
    )

dataset_info = info["dft_2d"]

print("Download URL:", dataset_info[0])
print("Internal JSON:", dataset_info[1])
print("Description:", dataset_info[2])
print("Reference:", dataset_info[3])
```

将结果写入：

```text
reports/database_source.txt
```

至少记录：

- 执行日期
- Python version
- jarvis-tools version
- 数据集名称
- 实际下载 URL
- ZIP 内部 JSON 文件名
- 数据描述
- 数据引用信息

原则：

**实际安装的软件和官方实时配置优先，不相信提示词中写死的旧 URL。**

---

# 10. Step 2：下载前验证 URL

从：

```python
get_db_info()["dft_2d"][0]
```

动态获取数据 URL。

使用：

```python
import requests

response = requests.get(
    url,
    stream=True,
    allow_redirects=True,
    timeout=60
)

print("Status:", response.status_code)
print("Final URL:", response.url)
print("Content-Type:", response.headers.get("Content-Type"))
print("Content-Length:", response.headers.get("Content-Length"))
```

保存：

```text
reports/download_probe.json
```

必须记录：

- HTTP status
- original URL
- final redirected URL
- Content-Type
- Content-Length
- redirect history

如果出现：

`403`

`404`

`500`

或者返回的是 HTML 页面而不是数据文件，则：

**立即停止。**

禁止假装下载成功。

---

# 11. Step 3：正式下载 dft_2d

使用官方接口：

```python
from pathlib import Path
from jarvis.db.figshare import data

raw_dir = Path("data/raw/jarvis")
raw_dir.mkdir(parents=True, exist_ok=True)

records = data(
    dataset="dft_2d",
    store_dir=str(raw_dir)
)

print("Object type:", type(records))
print("Number of materials:", len(records))
```

基本验收标准：

- `records` 应为 list；
- 材料数量应大于 1000。

不要写：

```python
assert len(records) == 1109
```

因为数据库以后可能更新。

但是如果出现：

`N < 900`

则需要停止并调查：

- 是否下载了错误数据集；
- 数据是否损坏；
- JARVIS 是否发生重大版本更新。

---

# 12. Step 4：保存原始数据库快照

原始文件下载后不得修改。

计算：

- file size
- SHA256
- modification time

例如：

```bash
sha256sum data/raw/jarvis/* > data/raw/jarvis/SHA256SUMS.txt
```

同时把 Python 解析后的数据库保存：

```text
data/raw/jarvis/dft_2d_snapshot.json
```

这样以后即使 JARVIS 更新，当前论文仍然可以完全复现。

---

# 13. Step 5：Schema Audit

不要根据论文猜字段。

首先获取所有真实字段：

```python
all_keys = sorted(
    set().union(*(record.keys() for record in records))
)

for key in all_keys:
    print(key)
```

针对每个字段统计：

- field
- dtype
- N_total
- N_nonmissing
- N_missing
- coverage_fraction
- N_unique
- sample_value

保存：

```text
data/audit/schema.csv
```

---

# 14. 缺失值定义

至少检查：

- None
- NaN
- ""
- "na"
- "NA"
- "N/A"
- "None"
- "null"
- "not available"

但是：

**禁止自动把 0 当作缺失值。**

例如：

`band_gap = 0`

很可能代表金属材料。

因此必须针对不同物性分别判断。

---

# 15. Step 6：自动搜索热电相关字段

对数据库所有字段名搜索以下关键词：

```text
seebeck
power
powerfact
conduct
conductivity
sigma
thermal
kappa
zt
mass
effective
elastic
gap
bandgap
boltz
carrier
mobility
exfol
ehull
formation
```

输出：

```text
reports/te_candidate_fields.txt
```

以下字段仅作为重点检查候选，不允许假定一定存在：

```text
jid
atoms
formation_energy_peratom
optb88vdw_bandgap
mbj_bandgap
elastic_tensor
effective_masses_300K
n-Seebeck
p-Seebeck
n-powerfact
p-powerfact
exfoliation_energy
ehull
```

实际数据库字段优先。

---

# 16. Step 7：建立 Property Coverage Matrix

输出：

```text
data/audit/property_coverage.csv
```

最终必须形成类似：

| Property | N total | N available | Coverage | Data form | Usable |
|---|---:|---:|---:|---|---|
| Structure | ? | ? | ? | Structure | ? |
| Band gap | ? | ? | ? | Scalar | ? |
| Effective mass | ? | ? | ? | Scalar/Vector | ? |
| n-Seebeck | ? | ? | ? | Scalar/Vector | ? |
| p-Seebeck | ? | ? | ? | Scalar/Vector | ? |
| n-PF | ? | ? | ? | Scalar/Vector | ? |
| p-PF | ? | ? | ? | Scalar/Vector | ? |
| Conductivity | ? | ? | ? | ? | ? |
| κ_e | ? | ? | ? | ? | ? |
| κ_L | ? | ? | ? | ? | ? |
| ZT | ? | ? | ? | ? | ? |
| Elastic data | ? | ? | ? | ? | ? |
| E_hull | ? | ? | ? | Scalar | ? |

所有问号必须通过实际数据库统计得到。

禁止人为猜测。

---

# 17. Property View 使用规则

根据真实数据量决定是否建立 property layer。

## A 类：主要物性层

满足例如：

`Coverage ≥ 50%`

并且：

- 单位明确；
- 计算条件明确；
- 数据物理意义明确。

则进入主要流形研究。

## B 类：Partial Property Layer

如果：

`N_available ≥ 100`

但覆盖率较低，

仍然建立 partial property graph。

## C 类：Exploratory Layer

如果：

`N_available < 100`

则只进行探索性分析。

100 和 50% 只是初始工作阈值。

后续必须进行 sensitivity analysis。

---

# 18. Step 8：确认 JARVIS 热电数据定义

对于数据库中实际存在的：

- n-Seebeck
- p-Seebeck
- n-power factor
- p-power factor
- conductivity
- κ_e
- κ_L
- ZT

必须通过 JARVIS 官方论文或官方 documentation 核实：

- temperature
- carrier concentration
- carrier type
- tensor direction
- relaxation time
- units
- Boltzmann transport method
- 2D normalization convention

写入：

```text
reports/thermoelectric_metadata.md
```

无法确定的项目写：

`UNVERIFIED`

禁止猜测。

---

# 19. 禁止自行反推数据库不存在的物性

理论关系：

`PF = S²σ`

但是除非可以确认：

- S 与 PF 来自相同温度；
- 相同载流子浓度；
- 相同载流子类型；
- 相同晶向；
- 相同 relaxation-time convention；
- 相同单位体系；

否则禁止使用：

`σ = PF / S²`

自行生成 conductivity。

默认设置：

```text
DERIVE_SIGMA = False
```

同样禁止自行通过：

`κ_e = LσT`

补充数据库不存在的电子热导率。

核心原则：

**只使用 JARVIS 中真正已经存在的数据。**

---

# 20. Step 9：验证所有二维结构

对于每个材料：

```python
from jarvis.core.atoms import Atoms

atoms = Atoms.from_dict(record["atoms"])
```

检查：

- 是否可以成功读取；
- n_atoms；
- elements；
- lattice；
- coordinates；
- vacuum direction。

输出：

```text
data/audit/structure_validation.csv
```

---

# 21. Step 10：建立 Material Master Table

创建：

```text
data/processed/material_master.parquet
```

使用：

`jid`

作为材料唯一身份 ID。

不要根据 chemical formula 去重。

例如：

`MoS2`

可以存在不同结构和不同 JID。

至少保留：

- jid
- formula
- elements
- n_atoms
- structure_valid
- number_of_available_views
- 数据库中真实存在的 scalar properties

---

# 22. 整体多流形框架

最终允许存在：

`M_struct`

`M_Eg`

`M_mstar`

`M_S`

`M_PF`

`M_sigma`

`M_kappa_e`

`M_kappa_L`

但：

**只有数据库真正存在相应数据时才建立。**

如果：

`N_available = 0`

则标记：

`SKIPPED_NOT_AVAILABLE`

不要为了匹配预设概念图而人为制造一个流形。

---

# 23. Step 11：建立结构流形

结构图记为：

`G_struct`

它作为整个多层材料网络的 backbone。

结构 descriptor 可以包括两部分。

## Composition descriptors

例如：

- 元素比例
- mean atomic mass
- atomic mass variance
- electronegativity mean
- electronegativity variance
- atomic radius statistics
- valence-electron statistics

## Structural descriptors

首选：

`SOAP`

同时测试：

- CrystalNNFingerprint
- coordination
- bond-length statistics
- bond-angle statistics
- symmetry
- layer thickness
- buckling
- area per atom

---

# 24. 二维材料真空层特殊处理

人为设置的真空层不能成为结构距离的重要来源。

例如同一个二维层：

`c = 15 Å`

和：

`c = 25 Å`

在结构空间中应该高度相似。

必须人为构造：

- 15 Å vacuum
- 20 Å vacuum
- 25 Å vacuum

三个测试结构。

如果 descriptor 明显改变，说明 descriptor 受到人工真空影响，需要修改结构表示方法。

---

# 25. 结构层禁止加入电子和热电性质

禁止将以下信息加入 `X_struct`：

- Seebeck
- PF
- band gap
- effective mass
- conductivity

`X_struct` 只允许包含：

**Composition + Geometry + Local Environment + Symmetry**

这样以后才能真正研究：

> 结构空间能否解释物性空间。

---

# 26. Step 12：建立 Property Layers

对于每一个实际存在的性质 v：

建立材料集合：

`D_v`

然后只在 `D_v` 中建立：

`G_v`

因此不同 property layer 的材料数量完全可以不同。

例如：

`N_struct = 1109`

`N_S = 900`

`N_PF = 850`

`N_mstar = 600`

都是允许的。

无需补全。

---

# 27. 一个标量不应过度称作高维流形

如果每个材料只有：

`S_i`

一个数，

那么这一层本质上属于一维性质相似关系。

此时更准确地称为：

**Seebeck similarity layer**

而不是宣称发现了复杂高维 Seebeck manifold。

如果每个材料拥有：

`X_S(i) = [S_n, S_p, S_xx, S_yy, S(T1), S(T2), ...]`

这样的多维响应向量，

才更适合称为：

`M_S`

即 Seebeck response manifold。

---

# 28. n 型和 p 型禁止平均

不要使用：

`S_avg = (S_n + S_p) / 2`

因为 Seebeck 可能一正一负，平均会造成没有物理意义的抵消。

优先分别建立：

`G_S_n`

和：

`G_S_p`

也可以测试联合表示：

`X_S = [S_n, S_p]`

两种方式都保留用于方法比较。

PF 同理：

`PF_n`

和：

`PF_p`

必须分别保存。

---

# 29. Step 13：每个 View 单独标准化

对于不同 `X_v`：

使用自己的 scaler。

禁止：

所有物性先拼接，再统一 StandardScaler。

普通连续变量可以使用：

`StandardScaler`

严重偏态数据可以测试：

`RobustScaler`

对于严格为正且跨多个数量级的数据，如：

- PF
- σ
- κ

优先测试：

`x_new = log10(x + epsilon)`

然后再标准化。

记录 epsilon 和全部 scaler 参数。

---

# 30. Step 14：建立 kNN Similarity Graph

对每个 view：

仅使用该 view 中具有真实数据的材料。

定义材料 i 和 j 在 view v 中的距离：

`d_v(i,j)`

建立 k-nearest-neighbor graph。

测试：

`k = 5, 10, 15, 20, 30, 50`

推荐局域尺度 Gaussian affinity：

`W_v(i,j) = exp[-d_v(i,j)² / (sigma_i × sigma_j)]`

其中：

`sigma_i`

可以取材料 i 的第 k 个近邻距离。

保存：

```text
graphs/W_structure.npz
graphs/W_bandgap.npz
graphs/W_mass.npz
graphs/W_Seebeck_n.npz
graphs/W_Seebeck_p.npz
graphs/W_PF_n.npz
graphs/W_PF_p.npz
```

只有真正存在的数据层才生成。

---

# 31. 每个 Layer 单独归一化

由于不同 property layers：

- 节点数不同；
- 图密度不同；
- 数值尺度不同；

禁止直接：

`W_total = W_1 + W_2 + W_3 + ...`

先分别进行：

- symmetric normalization

或：

- row-stochastic normalization

并明确记录方法。

---

# 32. Step 15：每个 View 独立 Embedding

每个足够大的 view 至少计算三类 embedding。

## PCA

作为线性 baseline。

## UMAP

主要用于二维可视化。

## Spectral Embedding 或 Diffusion Map

作为主要几何分析工具。

保存：

```text
manifolds/structure.csv
manifolds/seebeck_n.csv
manifolds/seebeck_p.csv
manifolds/pf_n.csv
manifolds/pf_p.csv
...
```

---

# 33. 禁止直接比较两个 UMAP 坐标

例如：

`UMAP_S(i) = (1.2, 3.1)`

和：

`UMAP_PF(i) = (1.2, 3.1)`

不能说明它们在两个流形中位置一样。

原因是不同 UMAP 的：

- origin 任意；
- rotation 任意；
- reflection 任意；
- scale 任意；
- global distance 可能严重失真。

因此不同 property views 必须通过：

**JID correspondence**

进行真正的 manifold alignment。

---

# 34. Step 16：建立 Partial Multilayer Graph

对于材料 i：

首先创建结构节点：

`M_i(struct)`

如果该材料有 Seebeck：

创建：

`M_i(S)`

如果有 PF：

创建：

`M_i(PF)`

如果有 effective mass：

创建：

`M_i(m*)`

如果缺少某个 property：

**不创建该节点。**

---

# 35. 使用 JID 作为跨层身份锚点

如果：

`M_i(struct)`

和：

`M_i(S)`

具有相同 JID，

则建立 identity edge：

`M_i(struct) ↔ M_i(S)`

同理：

`M_i(struct) ↔ M_i(PF)`

`M_i(struct) ↔ M_i(m*)`

等等。

这类跨层边表示：

**它们是同一个材料在不同物理观察空间中的表示。**

它不表示两个性质的数值相似。

---

# 36. 第一版采用 Star-Anchor 架构

推荐结构：

```text
                 Seebeck
                    |
                    |
Power factor --- Structure --- Band gap
                    |
                    |
              Effective mass
```

即：

**Structure layer = universal backbone**

暂时不直接建立：

`Seebeck ↔ PF`

`PF ↔ m*`

等所有两两跨层连接。

原因是：

- structure 覆盖率最高；
- 模型最清晰；
- 减少重复 identity constraints；
- 方便解释。

---

# 37. Supra-Adjacency Matrix

建立一个大矩阵 A。

其对角块是不同 layer 内部的材料相似图：

`W_struct`

`W_S`

`W_PF`

`W_mstar`

等等。

非对角块是相同 JID 的 identity correspondence。

可以概念性写成：

```text
A =

[ alpha_s * W_struct      lambda * C_sS      lambda * C_sPF      ... ]
[ lambda * C_sS^T         alpha_S * W_S      0                   ... ]
[ lambda * C_sPF^T        0                  alpha_PF * W_PF     ... ]
[ ...                     ...                ...                 ... ]
```

其中：

`alpha_v`

控制每个物理层内部关系的权重。

`lambda`

控制同一个材料在不同物理视图之间的身份约束强度。

---

# 38. 不要直接人为确定 lambda

第一轮可以：

`alpha_v = 1`

作为 baseline。

然后扫描：

`lambda = 0.01, 0.03, 0.1, 0.3, 1, 3, 10`

解释：

- lambda 很小时，不同 view 几乎独立；
- lambda 很大时，同一个材料的不同副本被强制重合。

寻找结果具有稳定性的 lambda 区间。

不要仅选择一组参数得到最好看的二维图。

---

# 39. Step 17：Joint Spectral Embedding

Supra adjacency matrix 为：

`A`

首先计算 degree：

`D_ii = Σ_j A_ij`

然后计算 normalized Laplacian：

`L_sym = I - D^(-1/2) A D^(-1/2)`

进行特征值分解。

去除 trivial eigenvector。

保留前 10–20 个 non-trivial eigenvectors：

`Phi = [phi_1, phi_2, ..., phi_d]`

作为统一多物理材料空间。

---

# 40. 公共空间不一定是二维

不要假设 intrinsic dimension = 2。

使用：

- eigenvalue spectrum
- eigengap
- participation ratio
- intrinsic-dimension estimation

判断合理的 latent dimension。

可能最终：

`d = 5`

或：

`d = 8`

才足够描述材料空间。

二维：

`(phi_1, phi_2)`

只用于画图。

真正进行：

- neighbor analysis
- distance calculation
- clustering
- cross-view tension
- material pair ranking

时，使用完整的 d 维 latent coordinates。

---

# 41. 最终每个材料只显示一个点

虽然 supra graph 中一个 JID 可能对应：

- structure node
- Seebeck node
- PF node
- effective-mass node

但最终统一材料图默认使用：

`M_i(struct)`

的 joint embedding 坐标。

这样保持：

**one JID = one point**

由于 structure node 已经通过 identity anchor 受到其他 property layers 的影响，

因此它已经不是纯结构坐标，而是：

**thermoelectric-informed material coordinate**

---

# 42. 同时保存所有 Property Node

不要删除其他 layer node 的 embedding。

保存：

```text
manifolds/all_multilayer_nodes.parquet
```

字段至少包括：

```text
node_id
jid
view
phi_1
phi_2
phi_3
...
```

这样可以研究同一个材料不同物理 representation 之间的偏移。

---

# 43. 定义 Cross-View Tension

对于材料 i 和 property view v：

定义：

`T_i,v = distance[Phi(M_i(struct)), Phi(M_i(v))]`

使用 d 维 Euclidean distance 或 diffusion distance。

如果：

`T_i,v ≈ 0`

说明材料的结构邻域和该物性邻域高度一致。

如果：

`T_i,v` 很大，

说明该材料属于：

**structure-property anomalous material**

即：

> 它的结构很像某类材料，但对应物性明显不像这一类材料。

这是非常值得关注的异常点。

---

# 44. 定义 Neighbor Overlap

对于两个物理 view：

`a`

和：

`b`

材料 i 各有 k 个近邻：

`NN_a(i)`

和：

`NN_b(i)`

定义：

`O_i(a,b) = |NN_a(i) ∩ NN_b(i)| / k`

然后对同时存在于两个 view 中的材料取平均：

`O_mean(a,b) = mean[O_i(a,b)]`

建立：

```text
              Struct   Eg   m*   S-n   S-p   PF-n   PF-p

Struct

Eg

m*

S-n

S-p

PF-n

PF-p
```

这样的 view similarity matrix。

---

# 45. View Similarity Matrix 的物理意义

它可以直接回答：

> 哪两个物理空间的材料邻域最一致？

例如如果：

`O_mean(struct, S)`

明显高，

说明：

**结构相似往往意味着 Seebeck 相似。**

如果：

`O_mean(struct, PF)`

很低，

说明：

**功率因子不能简单由结构邻域解释。**

这本身就是一个值得发表的材料物理结果。

---

# 46. 寻找 Structure-Close / Transport-Far 材料对

对于任意材料 i 和 j：

计算：

`d_struct(i,j)`

以及：

`d_transport(i,j)`

寻找满足：

`d_struct` 位于所有材料对最低 5%

同时：

`d_transport` 位于最高 20%

的材料对。

也就是：

**结构非常相似，但输运性质非常不同。**

保存：

```text
data/processed/structure_close_transport_far.csv
```

---

# 47. 为什么这些材料对对超晶格有价值

未来超晶格设计通常希望寻找：

**structural compatibility**

同时具有：

**electronic or transport contrast**

因此未来可以进一步寻找：

`d_struct ≪ 1`

`d_electronic ≪ 1`

但：

`d_thermal ≫ 1`

这意味着：

- 晶格结构相容；
- 电子性质相容；
- 热输运差异明显。

这类组合可能非常适合界面声子散射设计。

当前阶段只识别材料对，不生成超晶格。

---

# 48. Step 18：统一公共地图上的 Property Mapping

获得统一坐标：

`(Phi_1, Phi_2)`

后，

使用完全相同的坐标分别画：

- chemical family
- band gap
- effective mass
- n-Seebeck
- p-Seebeck
- n-PF
- p-PF
- E_hull
- exfoliation energy
- elastic properties

如果数据库实际存在：

- conductivity
- κ_e
- κ_L
- ZT

再加入。

不存在就跳过。

---

# 49. 缺失数据仍然显示在统一图中

例如绘制 Seebeck map 时：

有 Seebeck 的材料：

**按 Seebeck 数值着色。**

没有 Seebeck 的材料：

**灰色显示。**

不能直接删除这些材料。

这样同一张图可以同时表示：

- materials landscape
- property landscape
- property coverage

---

# 50. 定义材料信息完整度

定义：

`V_i = 材料 i 拥有的有效 property views 数量`

例如：

材料 A：

`V_A = 6`

材料 B：

`V_B = 3`

材料 C：

`V_C = 1`

可以：

- 点大小表示 V_i；
- 或单独绘制 information coverage map。

这样避免把只有结构信息的材料和具有大量输运数据的材料解释为相同置信度。

---

# 51. ZT 不能作为构建公共流形的输入

即使 JARVIS 中存在 ZT，

也禁止把它作为构建公共流形的输入特征。

因为：

`ZT = S²σT / (κ_e + κ_L)`

如果先把 ZT 输入流形，然后再声称发现 high-ZT region，

属于 target leakage。

ZT 只能作为：

**external validation label**

即：

公共流形建好以后再着色。

---

# 52. 如果没有可靠 κ_L

如果 Coverage Audit 发现：

`N_kappa_L = 0`

或数量极少，

则当前工作不要命名为：

**Full ZT Manifold**

而使用：

**2D Thermoelectric Transport Atlas**

重点研究：

`Structure → Electronic Structure → Seebeck → Power Factor`

之间的空间关系。

没有 κ_L 并不意味着这个工作无法进行。

---

# 53. Step 19：Property Smoothness

对于公共材料图中的相似度矩阵 W，以及材料性质 y：

定义：

`Smoothness(y) = Σ[W_ij × (y_i - y_j)²] / ΣW_ij`

分别计算：

- band gap
- effective mass
- Seebeck
- PF
- 其他实际存在的性质

然后随机打乱 property labels：

`1000 次`

得到 null distribution。

如果真实 Smoothness 显著低于随机情况，

说明：

**该物性确实沿材料流形连续变化。**

---

# 54. Step 20：流形稳定性验证

改变：

- kNN 的 k
- kernel bandwidth
- lambda
- latent dimension
- random seed
- SOAP parameters

重复计算。

至少测试约：

`20 组以上合理参数组合`

统计：

- trustworthiness
- kNN overlap
- Procrustes similarity
- graph-distance rank correlation

只有稳定存在的结构才能作为论文结论。

---

# 55. UMAP 只负责展示

UMAP 可以作为漂亮的二维可视化。

但以下分析禁止使用二维 UMAP Euclidean distance：

- material distance
- neighbor disagreement
- cross-view tension
- candidate pair ranking

这些分析主要使用：

- graph distance
- diffusion distance
- joint spectral coordinates

---

# 56. Random-Anchor Negative Control

保持每一个 property graph 内部结构不变。

然后随机打乱：

`JID correspondence`

重新进行 joint embedding。

如果随机 JID 对齐与真实 JID 对齐结果差别很小，

说明模型并没有真正利用跨物理空间对应关系。

因此 random-anchor test 必须作为负对照。

---

# 57. Baseline 1：Structure-Only

仅使用：

`G_struct`

建立：

`Phi_struct`

再与：

`Phi_joint`

比较。

分析：

> 加入电子和热电 property layers 后，哪些材料在公共空间中的位置变化最大？

这些材料可能具有强烈的 structure-property inconsistency。

---

# 58. Baseline 2：Complete-Case Concatenation

如果同时拥有所有主要 properties 的材料仍然足够多，例如：

`N_complete > 100`

则建立传统 baseline：

```text
complete samples
→ concatenate all features
→ scale
→ UMAP / spectral embedding
```

然后和 partial multi-view 方法比较。

目的：

证明 partial approach 可以使用更多材料和更多不完整信息。

---

# 59. Baseline 3：Independent Property UMAP

分别绘制：

- Structure UMAP
- Band-gap UMAP
- Effective-mass UMAP
- Seebeck UMAP
- PF UMAP

只作为：

**independent property landscapes**

不能直接把这些坐标拼接。

---

# 60. n 型和 p 型必须分别研究

分别绘制：

- S_n
- S_p
- PF_n
- PF_p

最终可以分别输出：

```text
top_n_candidates.csv
top_p_candidates.csv
```

不要默认最佳 n 型材料与最佳 p 型材料属于同一个材料区域。

---

# 61. 材料家族解释

对于公共材料空间不同区域，统计：

- elements
- chemical formula family
- symmetry
- layer group
- coordination
- band gap
- effective mass
- Seebeck
- PF

比较：

**Chemical Family**

**Structural Family**

**Transport Family**

三类材料分类是否一致。

如果出现明显不一致，

优先研究这些异常区域。

---

# 62. 重点筛选异常材料

对每个材料计算：

- cross-view tension
- neighbor disagreement
- local property residual
- number of available views

选取最异常的前 20–50 个材料。

人工检查：

- structure
- composition
- band gap
- effective mass
- Seebeck
- PF

这些异常点可能比单纯 PF 最大的材料具有更高科学价值。

---

# 63. 为未来超晶格研究输出材料 Pair Table

当前阶段禁止真正生成超晶格。

但可以输出：

```text
data/processed/candidate_pair_table.csv
```

至少包括：

```text
material_A
material_B

structure_distance
bandgap_distance
effective_mass_distance
Seebeck_distance
PF_distance

common_views
pair_information_completeness
```

未来再根据：

- lattice mismatch
- orientation
- band alignment
- phonon contrast
- interface direction

进一步研究 A/B 横向或垂直超晶格。

---

# 64. 当前阶段严格禁止

当前不进行：

- 异质结结构生成
- 超晶格结构生成
- VASP 新计算
- AMSET
- BoltzTraP
- phono3py
- ShengBTE
- ML 预测 κ_L
- ML 预测 mobility
- ML 缺失值填补
- 人工补充数据库

当前只利用：

**JARVIS 已经存在的数据。**

---

# 65. 建议论文 Figure 1：整体方法

```text
JARVIS dft_2d
        │
        ▼
Schema Audit
        │
        ▼
Property Coverage Audit
        │
        ▼
Partial Physical Views
        │
        ▼
View-Specific Similarity Graphs
        │
        ▼
JID Identity Anchors
        │
        ▼
Multilayer Supra Graph
        │
        ▼
Joint Spectral Embedding
        │
        ▼
2D Thermoelectric Materials Atlas
```

---

# 66. 建议 Figure 2：Structure-Only Map

首先建立纯结构空间。

分别按照：

- chemical family
- symmetry
- layer group

着色。

验证结构 descriptor 是否能够识别合理材料家族。

---

# 67. 建议 Figure 3：Independent Property Maps

分别展示：

- band-gap layer
- effective-mass layer
- Seebeck layer
- PF layer

比较不同物性空间中材料分类是否一致。

---

# 68. 建议 Figure 4：Unified Materials Atlas

使用完全相同的公共坐标：

`Phi_1`

和：

`Phi_2`

分别着色：

- E_g
- m*
- S_n
- S_p
- PF_n
- PF_p

这样不同 property landscape 可以直接视觉比较。

---

# 69. 建议 Figure 5：View Similarity Matrix

利用 Neighbor Overlap 建立矩阵。

回答：

> 哪些物理性质拥有最相似的材料邻域？

这是研究不同物理流形关系的核心图之一。

---

# 70. 建议 Figure 6：Cross-View Tension Map

将：

`T_i,v`

映射到公共材料图。

识别：

**structure-near / property-far**

异常材料。

---

# 71. 建议 Figure 7：超晶格前驱材料 Pair Map

横轴：

`d_structure`

纵轴：

`d_transport`

重点研究：

`d_structure 较小`

但：

`d_transport 较大`

的区域。

这些材料对可以作为后续超晶格研究的候选池。

---

# 72. 执行过程中必须设置 Checkpoint

禁止一次性生成整个项目并全部运行。

严格执行：

```text
Step 0
→ 编写程序
→ 执行
→ 检查 stdout/stderr
→ 检查输出文件
→ QA
→ 保存 checkpoint

Step 1
→ 编写程序
→ 执行
→ 检查
→ QA

Step 2
→ ...
```

发生异常时立即停止。

禁止：

```python
try:
    ...
except:
    pass
```

这种静默忽略错误的写法。

---

# 73. 当前第一轮只执行到 Coverage Audit

现在只执行以下阶段。

## Phase A

建立 Python 环境。

## Phase B

验证 `dft_2d` 当前真实官方下载来源。

## Phase C

真正下载 `dft_2d`。

## Phase D

保存原始快照并计算 SHA256。

## Phase E

读取全部数据库字段。

## Phase F

完成 Property Coverage Audit。

完成 Phase F 后：

**立即停止。**

当前不要开始流形分析。

---

# 74. 第一轮必须返回的结果

必须给出：

1. Python version

2. jarvis-tools version

3. `dft_2d` 实际下载 URL

4. 下载文件名

5. ZIP 内部 JSON 文件名

6. SHA256

7. 材料总数

8. 完整数据库字段列表

9. 每个字段的有效数据数量

10. 每个字段的 coverage

11. 热电相关字段及其数据结构

12. 每个字段的单位和计算条件是否能够确认

13. 哪些 property 可以进入下一阶段

最终形成：

| Property | N available | Coverage | Data form | Decision |
|---|---:|---:|---|---|
| Structure | ? | ? | structure | ? |
| Band gap | ? | ? | scalar | ? |
| Effective mass | ? | ? | scalar/vector | ? |
| n-Seebeck | ? | ? | scalar/vector | ? |
| p-Seebeck | ? | ? | scalar/vector | ? |
| n-PF | ? | ? | scalar/vector | ? |
| p-PF | ? | ? | scalar/vector | ? |
| Conductivity | ? | ? | ? | ? |
| κ_e | ? | ? | ? | ? |
| κ_L | ? | ? | ? | ? |
| ZT | ? | ? | ? | ? |

所有问号必须由真实下载后的数据库统计得到。

禁止猜测。

---

# 75. Coverage Audit 后再决定真正有哪些流形

不要预先假设一定存在：

`M_struct`

`M_S`

`M_sigma`

`M_kappa_e`

`M_kappa_L`

五个流形。

实际数据库可能最终只支持：

`M_struct`

`M_Eg`

`M_mstar`

`M_S`

`M_PF`

也可能还存在其他可用 property layers。

最终多流形体系必须：

**由真实数据库内容决定。**

---

# 76. 最终核心方法

整个研究可以概括为：

```text
Complete Structure Space
        +
Incomplete Electronic Property Spaces
        +
Incomplete Transport Property Spaces
        │
        ▼
JID-Based Partial Multi-View Alignment
        │
        ▼
Unified Materials Latent Space
```

也就是：

**一个完整的结构 backbone，加上若干不完整物理 property layers，通过相同 JID 作为身份锚点进行联合嵌入。**

---

# 77. 最终研究目标

这项研究不是：

> 从数据库中删除缺失行，再做一张 UMAP。

真正目标是：

> 最大程度利用一个天然不完备的二维材料数据库，分别学习不同物理性质下的材料邻域，并研究这些不同物理空间之间的几何关系。

重点回答：

**哪些结构变化会使二维材料从一个普通输运区域进入更好的热电输运区域？**

以及：

**哪些材料在结构空间中非常接近，但在输运空间中却显著不同？**

后一类材料尤其重要，因为它们将直接连接下一阶段的二维超晶格设计。

---

# 78. 当前立即执行的任务

现在从以下流程开始：

```text
Phase A
创建环境
    ↓
Phase B
验证 JARVIS dft_2d 数据源
    ↓
Phase C
下载数据库
    ↓
Phase D
验证文件完整性并保存快照
    ↓
Phase E
Schema Audit
    ↓
Phase F
Property Coverage Audit
    ↓
STOP
```

在真实的 Property Coverage Matrix 得到之前：

**不要进行任何流形建模。**

下一阶段的方法必须根据数据库中真正存在的 properties 和 coverage 决定。