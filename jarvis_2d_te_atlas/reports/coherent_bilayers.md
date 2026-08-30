# 共格双层超晶格构造与计算（Phase T 开端）

> 承接 42 的「晶格参数预筛配对」（89 对），本阶段真正构造共格双层超晶格并提交 MACE 计算稳定性与 κ_L。

## 一、构造方法（scripts/50_build_coherent_bilayers.py）

- 输入：JARVIS dft_2d 结构（修正的 `cart = frac @ lat` 行存储，`make_poscar` 同款）；配对于 `superlattice_candidate_pairs.csv`。
- 共格筛选：晶格失配 <1% 且 角度失配 <1° → 14 对。
- 构造：B 层旋转对齐 A 的平面内基矢 → 应变张量 `S = L_A @ inv(L_B)` 作用于 B 的平面内坐标（应变 = |S−I|fro）→ 层叠（层间范德华间隙 3.3 Å + 真空 15 Å）→ 写 POSCAR（ASE vasp5）。
- 过滤：平面内对齐（真空沿 z）、应变 <2%、层间距 3.0–4.6 Å。

## 二、构造结果（13 对成功，1 对 AsI3/Ge3O9Sb2 应变 2.2% 跳过）

| A (n) | B (p) | 应变 | 层间距 Å | 原子数 | c Å |
|---|---|---:|---:|---:|---:|
| OSn (JVASP-5888) | BrOSb (JVASP-77690) | 0.0001 | 3.30 | 10 | 30.2 |
| FHg (JVASP-28153) | BrOY (JVASP-76634) | 0.0004 | 3.82 | 10 | 40.6 |
| FHg (JVASP-28153) | BrOY (JVASP-60579) | 0.0006 | 3.82 | 10 | 40.6 |
| IOY (JVASP-60553) | BiBrO (JVASP-6217) | 0.0005 | 3.37 | 12 | 30.1 |
| OSn (JVASP-5888) | FHg (JVASP-28153) | 0.0047 | 3.81 | 8 | 31.1 |
| HoIO (JVASP-28173) | O3PbTi (JVASP-6244) | 0.0050 | 3.84 | 11 | 38.0 |
| BiBrO (JVASP-6217) | BrOSm (JVASP-20033) | 0.0052 | 4.32 | 12 | 31.3 |
| O3PbTi (JVASP-6244) | IOTm (JVASP-6259) | 0.0088 | 3.58 | 11 | 28.7 |
| O3PbTi (JVASP-6244) | BrOSm (JVASP-20033) | 0.0097 | 3.67 | 11 | 28.4 |
| HoIO (JVASP-28173) | BiBrO (JVASP-6217) | 0.0099 | 3.37 | 12 | 40.1 |
| BrH2Tb (JVASP-28115) | Br2V (JVASP-13546) | 0.0118 | 4.06 | 11 | 38.4 |
| FHg (JVASP-28153) | IOTm (JVASP-6259) | 0.0124 | 4.56 | 10 | 41.0 |
| FOV (JVASP-75331) | P (JVASP-77702) | 0.0130 | 3.80 | 8 | 34.7 |

POSCAR 与元数据：`data/superlattices/*_bilayer.POSCAR`、`coherent_bilayers_meta.csv`。

## 三、提交 MACE 计算（kl-mace-cpu 链：S1 优化 → S2 力 → S3 FC → S4 κ）

选 3 对代表（共格度最高 + n-p 输运对比强）：

| 材料 | 构成 | 应变 | |S_n−S_p| μV/K | A_ZT_e/B_ZT_e | 状态 |
|---|---|---:|---:|---:|---|
| OSn_BrOSb | SnO + BrOSb | 0.01% | 652 | 3.97/4.75 | S1 提交（3760883） |
| FHg_BrOY | FHg + BrOY | 0.04% | 714 | 6.98/3.71 | S1 提交（3760884） |
| IOY_BiBrO | IOY + BiBrO | 0.05% | 678 | 4.16/5.29 | S1 提交（3760899） |

- **稳定性**：S1 弛豫收敛 + 能量 vs 两个单独材料的 MACE 能量 → 界面能/结合能（S1 完成后计算）。
- **κ_L**：S2–S4 声子计算（κ 对称性审计随附）。
