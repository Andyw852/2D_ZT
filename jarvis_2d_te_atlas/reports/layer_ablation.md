# Layer Ablation（Q）

依次删除 n/p property 图谱中的一个 view，比较 consensus kNN(15) overlap（越低 = 该 layer 影响越大，λ=0.3）。

| carrier | drop | consensus overlap |
|---|---|---:|
| n | Eg | 0.217 |
| n | n_transport | 0.392 |
| n | m_electron | 0.575 |
| p | Eg | 0.227 |
| p | p_transport | 0.406 |
| p | m_hole | 0.588 |

结论：**Eg 是最强的 consensus driver**（金属/半导体二分主导公共空间），其次是 transport，mass 影响最小。
