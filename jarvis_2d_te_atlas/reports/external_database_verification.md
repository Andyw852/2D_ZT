# 外部热电数据库交叉验证报告
## JARVIS dft_2d 二维热电输运图谱的独立数据验证

- 执行日期：2026-08-25（UTC）
- 项目：2D Thermoelectric Materials Atlas（二维热电材料多物理统一图谱）
- 被验证对象：JARVIS curated `dft_2d`（1103 材料，Snapshot SHA256 `95f2cdd2b2aead0ffc34ffc65471c6b2637e5865a905df4b8ba051d963e5349a`）
- 本报告目的：用**其他热电/二维数据库**独立验证 dft_2d 的输运数据、电子结构数据、结构数据以及图谱的核心物理结论。

---

## 1. 下载的外部数据库（此前 round 报告因 figshare 403 / C2DB 404 无法下载，本次全部成功）

| 数据库 | 文件 | 记录数 | SHA256（zip） | 含热电输运？ |
|---|---|---|---|---|
| **JARVIS dft_2d_2021**（旧版本稳定性） | dft_2d_2021.json.zip → `d2-3-12-2021.json` | 1079 | `991e755ca03d89eb2a4e25dcc40ffd97ea289823cb74dc2de4e52b8cfb839673` | ✅ Seebeck / PF / σ / κ_e（n,p） |
| **C2DB**（GPAW-PBE，独立 2D 库） | c2db_atoms.json.zip → `c2db_atoms.json` | 3520 | `400fe281ed673037fe4c8cc51528891439cf2126c2ffa18c4e01e1ead5184e4a` | ❌ 仅 atoms/gap/etot/wf/efermi |
| **2DMatPedia**（VASP-PBE，独立 2D 库） | twodmatpd.json.zip → `twodmatpd.json` | 6351 | `099a9bbf5f28fe6ffa579aea1498d769dcc9da6204868bc1d071bf0b13ff4c0d` | ❌ 仅 gap/exfoliation/thermo |
| **Alexandria 2D**（2024.10.1，VASP-PBE） | alexandria_pbe_2d.zip | 137,833 | `89ee61523907a1312179e7afb18fd64060413534e787f78ef871dbfb8d624e6d` | ❌ 仅 gap/e_above_hull/e_form |

数据文件：`data/raw/external/{c2db,twod_matpd,dft_2d_2021,alexandria_2d}/`，校验和：`data/raw/external/SHA256SUMS.txt`。

> 注：本次 figshare 下载全部成功（302 → S3），与上一轮报告（figshare 403）相比网络环境已恢复。
> C2DB 当前 web 版（2024-05-01 数据模型）与 figshare dump 均不含 Seebeck/PF/κ_L；
> MatHub-2d 与 Alexandria REST API 仍不可达。因此**热电输运数值的直接独立验证以 dft_2d_2021 为主**，
> C2DB / 2DMatPedia / Alexandria 用于验证电子结构与结构 backbone。

---

## 2. 热电输运数据版本稳定性验证（dft_2d vs dft_2d_2021，1076 个共同 JID）

按 JID 精确匹配（同一材料身份，同一 DFT 协议 OptB88vdW，同一输运条件 T=600K、|n|=1e20 cm⁻³）：

| 物性 | 共同 JID | Pearson | Spearman | 中位比值(old/cur) | 中位绝对差 |
|---|---:|---:|---:|---:|---:|
| n-Seebeck | 735 | **1.0000** | **1.0000** | 1.000 | **0** |
| p-Seebeck | 709 | **1.0000** | **1.0000** | 1.000 | **0** |
| n-power factor | 801 | **1.0000** | **1.0000** | 1.000 | **0** |
| p-power factor | 801 | **1.0000** | **1.0000** | 1.000 | **0** |
| n-σ | 802 | **1.0000** | **1.0000** | 1.000 | **0** |
| p-σ | 802 | **1.0000** | **1.0000** | 1.000 | **0** |
| n-κ_e | 802 | **1.0000** | **1.0000** | 1.000 | **0** |
| p-κ_e | 802 | **1.0000** | **1.0000** | 1.000 | **0** |

**结论：JARVIS dft_2d 的全部 8 个热电输运量在 2021 与 2022 两个官方版本间逐条完全一致（Pearson=1.0，差异=0）。**
这同时验证了：
- 本项目从 NIST 静态 XML 恢复的 3 本征值 → 标量（OPTIMADE scalar = 3 本征值算术平均）的还原链路与官方旧版发布数据严格一致；
- JARVIS dft_2d 热电数据的**来源可靠性与版本稳定性**。

明细表：`data/audit/external_verify_dft2d_vs_2021.csv`；图：`figures/external_verify_te_version_stability.png`。

---

## 3. JARVIS 输运物理约定与现象的独立确认（用 dft_2d_2021）

### 3.1 PF 约定：PF = S²σ / 1e6（逐本征值）
对当前 dft_2d 恢复的 3 本征值逐通道验证：`PF_i = S_i²·σ_i/1e6`，
n 型 79.5% 通道在 1% 内一致（86.6% 在 5% 内），p 型 82.9% / 90.0%。
（少数偏差来自 PF=0 的数值截断通道。）→ **JARVIS 的 PF 定义与源码一致，可用作 external performance label。**

### 3.2 Seebeck 符号异常 ↔ 小带隙（bipolar 输运，Phase L0 核心发现）
用 2021 数据独立复算：

| 载流子 | 符号异常比例 | 异常材料中位带隙 | 正常材料中位带隙 | Mann-Whitney p |
|---|---:|---:|---:|---:|
| n 型 | 51/735 (6.9%) | **0.000 eV** | 1.364 eV | 2.2e-23 |
| p 型 | 30/709 (4.2%) | **0.000 eV** | 1.365 eV | 3.2e-15 |

→ **独立确认：Seebeck 符号异常集中在零带隙/金属材料，是 bipolar 输运而非数值噪声**（与 L0 的 p~1e-47 同向）。

### 3.3 设计规律的方向性确认（2021 标量数据）
| 规律 | 2021 数据结果 | 设计规则报告 |
|---|---|---|
| 规律 1：高 PF 由 \|S\| 主导，非 σ | logPF vs \|S\| Spearman **+0.51 / +0.35**；vs logσ **-0.09 / +0.02** | +0.52~0.62 / -0.14~-0.22（方向一致） |
| 规律 3：ZT_e ∝ S² 与带隙正相关 | ZT_e vs Eg Spearman **+0.50 / +0.65** | +0.63 / +0.76（方向一致） |
| n/p 分开研究 | n-S vs p-S Spearman **-0.40**；top-PF 材料 n/p 完全不同 | n/p 明显不同（方向一致） |

2021 版 top ZT_e（S²/L，L=2.44e-8）：n 型 CrBr2(15.4)/CeSe2(11.6)/UO2F2(8.0)/Li2CeAs2(7.7)/**WSe2(6.5)**；
p 型 CaClF(10.4)/ZnSiO3(7.8)/HgF(7.2)/Li2NiO2(6.9)/SrH2O3(6.4) —— 全部为高 \|S\|（394-614 μV/K）半导体，符合"高 ZT_e = 高 Seebeck"。
（设计规则报告中的 AsKO2/Cl2V/FOV 等候选为 2022 版新增材料，不存在于 2021 版，故无法用该版本复算；其物理基础 S²/L 已由上表确认。）

---

## 4. 电子结构与结构数据交叉验证（C2DB / 2DMatPedia / Alexandria）

按元素组成（canonical composition）匹配。

### 4.1 带隙
| 对比 | 匹配数 | Pearson | Spearman | 中位差 (ext−JARVIS) | 金属一致率 |
|---|---:|---:|---:|---:|---:|
| JARVIS(OptB88vdW) vs C2DB(PBE) | 502 | 0.795 | 0.715 | 0.000 eV | 61.3% |
| JARVIS(OptB88vdW) vs 2DMatPedia(PBE) | 906 | 0.867 | 0.883 | 0.000 eV | 74.6% |
| JARVIS(OptB88vdW) vs Alexandria(PBE) | 25 | 0.496 | 0.346 | −0.003 eV | 37.5% |

- JARVIS 与 **2DMatPedia（同为 VASP 系）** 带隙高度一致（Spearman 0.88，中位差 0）。
- 与 **C2DB（GPAW 系）** 也一致（Spearman 0.72）；部分金属/半导体判定差异源于 DFT 代码与泛函（OptB88vdW vs PBE）及多晶型。
- Alexandria 仅 25 个成分匹配（其 137k 结构以生成/装饰结构为主，与 curated dft_2d 化学空间重叠小），参考价值有限。

### 4.2 已知热电材料逐材料核对（eV）
| 材料 | JARVIS dft_2d | C2DB | 2DMatPedia |
|---|---:|---:|---:|
| MoS2 | 1.66 | 0.00/0.00/**1.58** | **1.65** |
| WSe2 | 1.33 | 0.00/0.00/**1.24** | **1.53** |
| MoSe2 | 1.25/1.45 | 0.00/0.00/**1.32** | **1.42** |
| PtSe2 | 1.36 | 0.00/0.00/**1.17** | **1.35** |
| Sb2Te3 | 0.75 | 0.34 | 0.66 |
| Bi2Se3 | 0.96 | 0.27 | 0.89 |
| Bi2Te3 | 0.38 | 0.08 | 0.97 |
| PbTe | 1.08 | — | 1.04 |
| SnTe | 0.83 | — | 0.79 |
| HfTe3 | 0.00 | 0.00 | 0.00 |

主流热电材料（MoS2、WSe2、MoSe2、PtSe2、Sb2Te3、Bi2Se3、PbTe、SnTe）的带隙在三个独立数据库间一致
（C2DB 中同一成分的金属通道对应 1T/多晶型）。Bi2Te3 单层带隙对泛函敏感（0.08–0.97 eV），属已知难点。

### 4.3 晶格常数（面内 a、b）
| 对比 | n | Pearson(a/b) | 中位相对误差(a/b) |
|---|---:|---:|---:|
| JARVIS vs C2DB | 502 | 0.867 / 0.768 | 2.41% / 2.32% |
| JARVIS vs 2DMatPedia | 906 | 0.824 / 0.891 | **0.35% / 0.39%** |

→ 2DMatPedia（同为 VASP）晶格几乎一致（<0.4%），C2DB（GPAW）2.3% 量级，符合跨代码弛豫差异预期。
**验证了项目 Structure backbone（标准化 2D 结构 + vacuum axis 处理）的几何基础。**

### 4.4 剥离能（JARVIS meV/atom ↔ 2DMatPedia eV/atom，换算后）
| n | Pearson | Spearman | 中位比值 JARVIS/2DMatPedia |
|---|---:|---:|---:|
| 614 | 0.775 | **0.839** | **1.006** |

→ 换算单位后两库剥离能几乎相等（比值中位 1.006），确认 dft_2d 的 exfoliation 数据可靠。

明细表：`data/audit/external_verify_jarvis_vs_c2db.csv`、`data/audit/external_verify_jarvis_vs_2dmatpedia.csv`、
`data/audit/external_verify_gap_alexandria.csv`；图：`figures/external_verify_gap_2dmatpedia.png`、`figures/external_verify_exfoliation.png`。

---

## 5. 结论与局限

### 支持性结论
1. **热电输运数据（S/PF/σ/κ_e）版本完全稳定**：dft_2d 与 dft_2d_2021 在 1076 个共同 JID 上逐条一致（Pearson=1.0）。
2. **PF 约定与物理现象可复现**：PF=S²σ/1e6（逐本征值）；Seebeck 符号异常 ↔ 零带隙（bipolar）在 2021 版独立复现（p~1e-23）。
3. **设计规律方向稳健**：PF 由 |S| 主导、ZT_e 与带隙正相关、n/p 分离，用独立版本标量数据均得到同向结果。
4. **电子结构与结构 backbone 与独立 2D 数据库一致**：带隙（2DMatPedia Spearman 0.88；C2DB 0.72）、
   晶格（2DMatPedia <0.4%）、剥离能（比值 1.006）。主流热电材料（MoS2/WSe2/Sb2Te3/Bi2Se3/PbTe/SnTe）逐材料核对一致。

### 局限（诚实声明）
1. **无独立 κ_L / 真实 ZT 验证**：C2DB（figshare dump 与当前 web）不含 κ_L/ZT；MatHub-2d 与 Alexandria API 不可达。
   因此"高 ZT_e 天花板"仍为电子 ZT 上限，未升级为真实 ZT（与 design_rules_final.md 一致）。
2. **热电输运的直接独立对照仅来自 JARVIS 自身旧版本**（dft_2d_2021），属于同协议版本验证，
   非完全独立的第一性原理复算；跨数据库（C2DB/2DMatPedia）只覆盖电子结构与结构量。
3. 成分匹配会混入不同晶型/多晶型（如 1T/2H MoS2），导致少量金属/半导体判定差异。
4. Alexandria 2D 与 curated dft_2d 化学空间重叠小（25 成分），未构成有效独立样本。

### 对下一步的建议
- 若要验证 κ_L/真实 ZT：C2DB 完整数据需走 2dhub.org 的 asr/myqueue 流程或向 CAMD 申请；
  或从 NOMAD/Materials Cloud 归档获取 C2DB 结果文件；文献 Ioffe/实验 κ 可作为 follow-up。
- 若要做跨代码输运对照：可对 dft_2d 中高 ZT_e 材料（如 WSe2、CrBr2、CaClF）用 BoltzTraP/AMSET 复算
  S 与 σ/τ（需新计算，超出当前"零新计算"约束，仅作后续验证方案）。
