# 参考 autonomous-research-engine 的框架改进分析

> 被参考对象：`autonomous-reasearch-plugin-system.zip`（确定性多智能体研究引擎，DDM=design/compute/measure）。
> 改进对象：`matexplore/`（我的 4-agent 闭环：analyst→generator→validator→collector）。

## 一、ARPS 值得借鉴的 6 个点

1. **机制优先的假设语法**（P1）：每条假设必须写成 `设计变量 → 可算描述子链 → 目标性质 → 方向`，
   纯统计相关（如"ρ>0.7"）不允许作为 claim。→ 我现在的假设就是纯 Spearman 相关，需要升级。
2. **对抗性审稿（Challenger 角色）**：独立批评者攻假设的边界/混淆/替代机制/决策价值。
   → 我的框架完全没有这一步，是最缺的一块。
3. **对比先于相关**（P4）：单变量受控对比决定机制，描述子相关只"暗示"机制。
   → 我目前只有相关 + 原型替换，缺"seed vs 替换后"的受控对比设计。
4. **证据分级**（P5）：screening(粗筛)/discrimination(甄别)/confirmation(确认)，方法与决策等级匹配。
   → 我的验证层没标证据等级，容易把粗筛当确认。
5. **规律成熟度阶梯**（lead→rule_candidate→rule→law）+ 研究前沿（frontier），
   负结果/失败记录必须保留。
6. **职责分离**：假设作者≠证据审稿人；执行者≠结果解释者。

## 二、本轮落地的改进

- 新增 **Challenger 智能体**（对抗性审稿，7 问 + 推荐甄别性实验）；
- Analyst 假设改为**机制优先**（`变量→描述子→目标→方向` + 边界 + 可证伪条件 + next_best_test）；
- Validator 给每条结果打**证据等级**（screening/discrimination/confirmation）；
- Ledger 增加 **frontier/规律成熟度** 记录；
- 生成器的单元素替换被显式框为**单变量受控对比**（seed vs 替换后）。

## 三、关键发现（Challenger 会据此出问题）

- 组成特征（电负性/电离能）与 ZT_e 的 +0.34~+0.39 很可能**被 Eg 混淆**（Eg 才是主驱动 +0.63~0.76）；
- 我的单元素替换同时改变**电负性 + 原子质量 + 原子半径**（如 Pb→Sn），不是严格单变量 → 需明确"多描述子联动"；
- 甄别性实验应是：**对 seed 与替换后候选算 Eg（或 m\*）**，看方向是否与预测一致，而不是再看一个相关。

