# JARVIS dft_2d 热电数据定义核实（Step 8，已按源码修正）

> 原则：只写能确认的；无法确定的项目一律写 `UNVERIFIED`，禁止猜测。

## 核心结论（直接回答「scalar / tensor / 固定温度与 doping」）

dft_2d 的热电输运量 **S（Seebeck）、σ（电导率）、κ_e（电子热导）既不是单个标量，也不是 3 个原始 tensor 分量，
而是 3×3 张量对角化之后的 3 个本征值（eigenvalues），且来自单一固定温度、单一固定掺杂浓度。**

证据链（jarvis-tools 2026.6.12 源码 `jarvis/db/vasp_to_xml.py::boltztrap_data`）：

```python
def boltztrap_data(self, path=..., temperature=600, doping=1e20):
    ...
    small_p = all_data["condtens_fixdoping"]["p"][temperature]   # 固定 T = 600 K
    small_n = all_data["condtens_fixdoping"]["n"][temperature]
    for i, j in small_p.items():
        if j["N_cm3"] == doping:                                  # 固定 n = 1e20 cm^-3
            pseeb  = eigvals(tmp["seeb"].reshape(3,3) * 1e6)      # S: 3×3 对角化 → 3 本征值, ×1e6 → μV/K
            pcond  = eigvals(tmp["cond"].reshape(3,3)) / 1e14     # σ: 3 本征值 / 1e14（常数 τ 约定）
            ppf    = pseeb**2 * pcond / 1e6                       # PF = S²σ/1e6（逐本征值）
            pkappa = eigvals(tmp["kappa"].reshape(3,3))           # κ_e: 3 本征值（原始 κ_e/τ，未缩放）
```

调用点 `main_boltz_data()` 里 `self.boltztrap_data(path=path)` **未传 temperature/doping**，即使用上述默认值
**T = 600 K、doping = 1e20 cm^-3**。

### 三条子问题的答案

1. **是 scalar 吗？** 不是。原始 figshare `d2-12-12-2022.json` 里每个量是 **3 个本征值**（3 分量）。
   本工作使用的 OPTIMADE 导出把这 3 个本征值**进一步取平均成了单个标量**（所有值都带 ÷3 尾数，如
   `n-Seebeck = -101.86333333333333 = (-305.59)/3`，`nkappa×3 = 993086158000000.0` 恰为整数，实证为 3 值平均）。
2. **是 3 分量 tensor 吗？** 是「3 分量」，但不是原始的 xx/yy/zz 或任意坐标系分量，而是 **3×3 张量对角化后的 3 个本征值**
   （主轴方向的主值）。原始 9 个分量（含非对角）在 `boltztrap.condtens_fixdoping` 里，入库前被 `np.linalg.eigvals` 压缩为 3 个。
3. **是固定温度 / 固定 doping 吗？** 是。**单一 T = 600 K、单一 |n| = 1e20 cm^-3**（n 型取 -1e20，p 型取 +1e20），
   不是温度扫描、也不是 doping 扫描后平均。

## 逐字段核实表（按源码）

| 字段 | 物理量 | 温度 | 掺杂 | 载流子类型 | 数据结构（原始→本工作） | 弛豫时间 | 单位 | 结论 |
|---|---|---|---|---|---|---|---|---|
| `_jarvis_n-Seebeck` | n 型 Seebeck | **600 K（固定）** | **1e20 cm⁻³（固定）** | n | 3 本征值 → 均值(标量) | 常数 τ（S 与 τ 无关） | **μV/K（源码 ×1e6 确认）** | 可用 |
| `_jarvis_p-Seebeck` | p 型 Seebeck | 600 K | 1e20 | p | 3 本征值 → 均值 | 同上 | **μV/K 确认** | 可用 |
| `_jarvis_n-powerfact` | n 型功率因子 | 600 K | 1e20 | n | 3 本征值 → 均值 | 常数 τ | PF = S²σ/1e6（UNVERIFIED 最终单位） | 谨慎 |
| `_jarvis_p-powerfact` | p 型功率因子 | 600 K | 1e20 | p | 3 本征值 → 均值 | 常数 τ | 同上 | 谨慎 |
| `_jarvis_ncond`/`_jarvis_pcond` | 电导率 σ | 600 K | 1e20 | n/p | 3 本征值 → 均值 | **常数 τ（intrans tau=0 确认）** | 源码 `/1e14`（≈ τ=1e-14 s 约定下的 σ，UNVERIFIED） | 谨慎 |
| `_jarvis_nkappa`/`_jarvis_pkappa` | 电子热导 κ_e | 600 K | 1e20 | n/p | 3 本征值 → 均值 | 常数 τ | **原始 κ_e/τ（未缩放，故数值 ~1e13–1e14）** | 谨慎/暂缓 |
| `_jarvis_avg_elec_mass` | 电子有效质量 | — | — | — | 3 值均值 | — | UNVERIFIED（数值 0.003–59000，疑含 2D 归一化） | 谨慎 |
| `_jarvis_avg_hole_mass` | 空穴有效质量 | — | — | — | 3 值均值 | — | UNVERIFIED | 谨慎 |

## 关键结论

1. **Seebeck 单位确认为 μV/K（源码 ×1e6），覆盖率 73.1%/72.7%，可直接进入主要物性层。**
2. **所有输运量固定 T=600 K、|n|=1e20 cm⁻³、常数弛豫时间**。这是「一个条件点」的数据，
   不是温度/doping 响应面；因此每个材料每个量在原始库里只有 3 个数（3 个本征值），本工作的 OPTIMADE 标量是它们的均值。
3. **σ 与 κ_e 采用常数 τ**：σ 被源码 ÷1e14（等价于假定 τ=1e-14 s），κ_e 未缩放（原始 κ_e/τ）。
   绝对数值不能直接物理解释为材料本征 σ 或 κ_e；用于同一计算协议内的相似性/邻域/排序仍可用。
4. **禁止反推** σ = PF/S² 或 κ_e = LσT（τ 约定与缩放不一致）。
5. **κ_L 不存在（N=0）、ZT 不存在（N=0）** ⇒ 本工作命名为 **2D Thermoelectric Transport Atlas**，
   研究 `Structure → Electronic Structure → Seebeck → Power Factor`。
6. **弹性张量在 OPTIMADE 导出中全为 -99999/None（N=0）** ⇒ `SKIPPED_NOT_AVAILABLE`。
7. 数据源说明：figshare（jarvis-tools 默认）在本环境网络层 403；改用官方 NIST OPTIMADE API
   `https://jarvis.nist.gov/optimade/jarvisdft/v1/structures/`（entry id 前缀 `dft_2d_`）。
