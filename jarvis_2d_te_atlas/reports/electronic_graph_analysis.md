# Electronic Graph 分析（Phase N）

## Eg 层（G_Eg_v1, N=1103）

- 只使用 Eg_optb88vdw，距离 = |Eg_i - Eg_j|（Eg=0 是金属，不是缺失）。
- 397 金属（Eg=0）+ 706 半导体。
- 图在任意 k 下都稳定地分裂为 **2 个连通分量**：金属 clique（397）与半导体链（706），giant=0.64。
- 这是物理上正确的结果（band-gap similarity layer 天然分金属/半导体），不是图构建缺陷。

## Electronic-n/p 层（N=678）

- Electronic-n = [Eg, m_elec_median]，Electronic-p = [Eg, m_hole_median]，RobustScaler + Euclidean。
- 678 材料中 204 金属 + 474 半导体，同样分裂为金属/半导体两个分量，giant=0.70。
- Electronic-n 与 Electronic-p 距离相关 0.74（共享 Eg）。

## 关键结论

- effective mass 数据不完整（678/1103），未做填补；Eg 层（1103）与 rich 层（678）分层保留。
- 金属/半导体二分是 Electronic View 的本质结构，后续靠 Structure backbone（全连通）桥接。
