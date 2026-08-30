# JARVIS dft_2d 第一轮（Phase A–F）执行报告

执行日期：2026-08-24（UTC）
项目：2D Thermoelectric Materials Atlas（二维热电材料多物理统一图谱）

## 1. 环境

- Python：3.11.16（conda env `te_manifold`）
- jarvis-tools：2026.6.12
- 关键依赖：numpy 2.4.6 / pandas 3.0.5 / scipy 1.17.1 / scikit-learn 1.9.0 / pymatgen 2026.5.4 / dscribe 2.1.2 / umap-learn 0.5.12 / pyarrow 25.0.1 / networkx 3.6.1 / ase 3.29.0
- 完整见 `reports/pip_freeze.txt`

## 2. 数据源（重要）

- jarvis-tools 官方配置 `get_db_info()["dft_2d"]`：
  - URL：`https://ndownloader.figshare.com/files/38521268`
  - ZIP 内部 JSON：`d2-12-12-2022.json`
- **figshare 在本执行环境被网络层 403 阻断**（figshare.com / api.figshare.com / ndownloader.figshare.com 全部 403，
  见 `reports/download_probe.json`），因此**未假装下载成功**，改用官方 NIST 托管的
  **JARVIS-DFT OPTIMADE API**（`https://jarvis.nist.gov/optimade/jarvisdft/v1/structures/`）作为权威数据源，
  按 entry id 前缀 `dft_2d_`（`_jarvis_source == "dft_2d"`）精确筛选 2D 数据。

## 3. 快照与校验

- 快照文件：`data/raw/jarvis/dft_2d_snapshot.json`（3.964 MB）
- SHA256：`95f2cdd2b2aead0ffc34ffc65471c6b2637e5865a905df4b8ba051d963e5349a`
- 材料总数：**1103**（unique JID 1103，> 1000 验收通过）

## 4. Schema

- 字段总数：**78**（每条记录的 `attributes` 内；含标准 OPTIMADE 结构字段 + `_jarvis_*` 物性字段）
- 完整字段与逐字段 coverage 见 `data/audit/schema.csv`

## 5. Property Coverage Matrix（真实统计，非猜测）

| Property | N available | Coverage | Data form | Decision |
|---|---:|---:|---|---|
| Structure | 1103 | 100% | structure | A |
| Band gap (OptB88vdW) | 1103 | 100% | scalar (eV) | A |
| Band gap (MBJ) | 246 | 22.3% | scalar (eV) | B |
| Band gap (HSE06) | 54 | 4.9% | scalar (eV) | C |
| Effective mass (elec/hole) | 678 | 61.5% | scalar | A* |
| n-Seebeck | 806 | 73.1% | scalar (μV/K) | A |
| p-Seebeck | 802 | 72.7% | scalar (μV/K) | A |
| n-Power factor | 802 | 72.7% | scalar | A* |
| p-Power factor | 802 | 72.7% | scalar | A* |
| n/p-Conductivity | 802 | 72.7% | scalar | A* |
| κ_e (n/p, electronic) | 802 | 72.7% | scalar | A*（单位 UNVERIFIED） |
| Dielectric (static eps) | 887 | 80.4% | scalar×3 | A |
| Exfoliation energy | 748 | 67.8% | scalar | A |
| E_hull | 1103 | 100% | scalar | A（全 0.0） |
| Elastic tensor | 0 | 0% | tensor | SKIPPED |
| **κ_L (lattice)** | **0** | **0%** | — | **SKIPPED_NOT_AVAILABLE** |
| **ZT** | **0** | **0%** | — | **SKIPPED_NOT_AVAILABLE** |

* 采用常数弛豫时间，绝对数值尺度不可直接物理解释，用于相似性/邻域仍可用。

## 6. 下一阶段可进入的物性层

- 主要层（A）：Structure、OptB88vdW 带隙、有效质量、n-Seebeck、p-Seebeck、n-PF、p-PF、n/p 电导率、静态介电、剥离能、E_hull
- Partial 层（B）：MBJ 带隙（N=246）、高频介电（N=249）
- Exploratory（C）：HSE06 带隙（N=54）
- 跳过：Elastic tensor、κ_L、ZT

## 7. 结论

- 本数据源下 **κ_L 与 ZT 均不存在** ⇒ 本工作应命名为 **2D Thermoelectric Transport Atlas**，
  重点研究 `Structure → Electronic Structure → Seebeck → Power Factor` 的空间关系。
- n/p 型 Seebeck 与 PF 需分别研究（禁止平均）。
- 单位/条件核实详见 `reports/thermoelectric_metadata.md`（无法确认者已标 UNVERIFIED）。

—— 按第 73 节要求，Phase F 完成后 **STOP**，未进行任何流形建模。
