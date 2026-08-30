# Random Anchor Negative Control（Q）

保持所有单视图图不变，仅打乱 JID identity correspondence，重算 joint embedding，比较真实 anchor 与随机 anchor 的 transport preservation（λ=0.3）。

| carrier | P_transport real | P_transport null | z | p |
|---|---:|---:|---:|---:|
| n | 0.607 | 0.558±0.011 | 4.4 | 0.000 |
| p | 0.637 | 0.584±0.009 | 5.8 | 0.000 |

结论：**真实 JID anchor 显著优于随机 anchor**（z=4.4/5.8，p=0.000）。
- JOINT_ALIGNMENT_VALID = True。
- 注意：anchor tension（绝对）真实 0.0466 略高于随机 0.0406，但 preservation 才是判断"对齐是否有意义"的正确指标——真实 anchor 显著更好地保留了 view 内几何。
