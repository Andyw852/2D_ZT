# 不同特征流形之间的距离图

## 图的含义

这不是把不同 UMAP 图的二维坐标直接相减。每个特征视图先在同一批材料上产生完整的材料—材料距离矩阵，随后定义两个流形的距离为：

\[
d(A,B)=\sqrt{2\left[1-\rho_s\left(\mathrm{vec}(D_A),\mathrm{vec}(D_B)\right)\right]},
\]

其中 \(\rho_s\) 是两个流形内部材料对距离的 Spearman 相关。最后仅对视图之间的距离矩阵做 classical MDS。图中横纵轴没有单独物理含义，点云间距才有意义。

- 左图：Structure、Band gap、Electronic n/p、Transport n/p 在 674 个完整共有材料上的距离。
- 右图：Structure、Band gap、Electronic、Elastic、Lattice κL 在 58 个五视图完整材料上的距离。
- 小点/椭圆：300 次、每次抽取 80% 材料后重新计算的位置，表示样本不确定性。

## 主要结果

1. n/p 电子流形最接近：\(\rho=0.741\)，\(d=0.719\)。
2. n/p 输运流形也共享明显结构：\(\rho=0.443\)，\(d=1.055\)。
3. 结构与电子/输运流形均很远：\(\rho=0.034\)–0.083，说明结构近邻不能直接替代电子输运近邻。
4. κL 子集中，κL 与弹性（\(\rho=0.227\)）和结构（\(\rho=0.192\)）比与电子流形（\(\rho=0.040\)）更接近。
5. κL 子集只有 58 个完整材料，因此 Structure/Elastic 云团明显更宽；这应解释为结论不确定性，而不是额外的物理维度。

## 产物

- 静态图：`figures/feature_manifold_distance_map.png`、`figures/feature_manifold_distance_map.pdf`
- 数值表：`data/audit/manifold_of_manifolds_distances.csv`
- 可视化数据：`data/processed/manifold_of_manifolds_map.json`
- 复现脚本：`scripts/46_plot_manifold_of_manifolds.py`
