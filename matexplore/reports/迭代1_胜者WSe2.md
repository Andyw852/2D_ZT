# 迭代 1 结果：胜者 = WSe2

计算平台：带隙 unihamgnn(3090) + MACE 链(opt/phonon) 已切到 jzzn。

## 完整验证链结果（已跑通的材料）

| 材料 | 带隙 eV (unihamgnn SOC) | E/atom eV (MACE) | 声子稳定 | 判定 |
|---|---|---|---|---|
| **WSe2** | **1.166** | **-7.210** | **stable(0 虚频)** | ✅ **胜者** |
| InSe | 0.250 | 待 | 待 | 窄带隙，边缘 |
| Bi2Te3 | -0.757 (金属) | -3.458 | stable | ❌ 金属，非 TE 半导体 |
| Br2ClFPb2 | 2.576 (太宽绝缘体) | -3.439 | unstable(-4.2THz) | ❌ 不稳 + 绝缘体 |
| 其余 8 | 待 | 待 | 待 | 生成候选，多数已在此前 MACE 声子中被判不稳 |

## 为什么 WSe2 是胜者

1. **带隙 1.166 eV**：理想热电半导体窗口（约 0.3–1.5 eV），不是金属(Bi2Te3 -0.76)也不是绝缘体(Br2ClFPb2 2.58)。
2. **动力学稳定**：jzzn MACE+phonopy 声子谱 0 虚频（min_freq -0.0003 THz），是可真实存在的二维相。
3. **MACE 弛豫收敛**：E/atom -7.21 eV，力 4.6e-4 eV/Å、应力 0.001 GPa，处于势能面极小点。
4. **与 JARVIS 库一致**：WSe2 本就是 JARVIS 中 n 型高 ZT_e 的已知二维热电材料，验证了整条「生成→带隙→稳定性」流水线找回了物理上正确的答案。

## 数据出处

- 带隙：unihamgnn(Uni-HamGNN 通用 SOC 哈密顿量，3090)，band_summary.json。
- 稳定性：MACE-matpes-pbe-omat-ft + phonopy(3×3×1)，jzzn 登录节点。
- 结果文件：matexplore/reports/iteration1_bandgap.csv。

## 剩余步骤

- KL(晶格热导率)是真实 ZT 的最后一块，用 kl-mace-cpu(jzzn) / kl-mace-gpu(3090) 可算；
- 其余 8 个候选的带隙/声子需在 3090 恢复后重试 unihamgnn，或逐个在 jzzn 跑 MACE 链。
