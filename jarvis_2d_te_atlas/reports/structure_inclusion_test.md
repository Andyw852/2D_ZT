# Structure Inclusion Test（P-F）

| 配置 | λ=0.3 P_transport | λ=1.0 P_transport | P_structure |
|---|---:|---:|---:|
| n-Full（含 Structure） | 0.587 | 0.645 | 0.46 |
| n-Property（无 Structure） | 0.607 | 0.555 | — |
| p-Full | 0.616 | 0.664 | 0.46 |
| p-Property | 0.637 | 0.554 | — |

结论：
1. Structure inclusion 对 property preservation 影响**大致中性**（±0.05），不"严重破坏"。
2. 但 **Structure 自身 geometry 几乎无法保留**（P_structure 0.13–0.51，始终 < 0.60）。
3. 根因：Structure 与 property views 的 distance correlation ≈ 0.02–0.07（去相关）。
4. 决策：**Structure = linked auxiliary view，不进公共空间**（不是因为它破坏 property，而是因为它与 property 本质去相关、无法有意义对齐）。
