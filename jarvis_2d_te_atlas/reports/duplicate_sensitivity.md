# Duplicate Sensitivity（Q）

- 存在 d_structure = 0 的精确重复（不同 JID 同结构），共 **12 个 exact-duplicate 组**。
- 主模型保留全部 JID；**真正重跑** collapsed 联合流形（每组用代表 JID，共 1091 个代表）并比较 transport kNN(15) 邻域保持：
  - n 型：collapsed kNN overlap = **0.960**（N=796）
  - p 型：collapsed kNN overlap = **0.921**（N=793）
- **DUPLICATE_SENSITIVITY_PASSED = True**（去重后代表材料的 transport 邻域保持 0.92–0.96，重复条目未显著歪曲主要科学结论）。
