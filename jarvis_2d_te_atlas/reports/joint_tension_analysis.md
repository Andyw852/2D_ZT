# Joint Cross-View Tension（R）

T_i,v = distance(Consensus_i, ViewCopy_i,v)，在 λ=0.3 的 property 图谱 20 维 joint space 中计算。

| carrier | T_transport 中位 | T_Eg 中位 | T_mass 中位 |
|---|---:|---:|---:|
| n | 0.0418 | 0.0604 | — |
| p | 0.0422 | 0.0640 | — |

高 T_transport（cross-view transport anomaly，top-5）：
- n：JVASP-14456, JVASP-27864, JVASP-7033, JVASP-27853, JVASP-20012。
- p：JVASP-13619, JVASP-5947, JVASP-6601, JVASP-6805, JVASP-786。

（完整 top-50 见 data/processed/joint_tension_{n,p}.csv；高 T 只能称 cross-view anomaly，不称 superlattice candidate。）
