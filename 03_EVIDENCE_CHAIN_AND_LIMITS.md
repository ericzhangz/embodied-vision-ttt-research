# 证据链、版本归属与结论边界

> **历史证据链，正文保留。** 本文记录旧独立/共享门控阶段；当前源码已包含关联能量方案。新证据请读 [T1 实现记录](references/RELATION_ENERGY_IMPLEMENTATION_20260909.md)和[匹配修正检验](references/RELATION_ENERGY_SELECTIVITY_20260909.md)，不沿用本文的“源码仍是独立gate”作为当前状态。

## 1. 怎样定位与复算

源码和原始量优先于摘要、测试数和叙事。[源码快照](code/conductance_star_arithmetic_check.js)仍是独立gate版本；当前query共享gate尚未实现，包中没有aligned运行JSON。对保存向量做内积/逐坐标乘法属于解析复算，不是执行新模型。完整性以[MANIFEST](MANIFEST.json)为准。

下文指针使用JSON Pointer：数字是从0开始的数组下标。为避免反复打印长路径，先定义精确前缀，后缀直接拼接；不是模糊搜索关键词。

| 代号 | 文件及根指针 |
|---|---|
| R | [互易电流算术](evidence/CONSTRUCTIVE_COUPLING_CHECK_20260909.json)，根 `/constructiveCouplingQualification` |
| S | [三路后继关联](evidence/STRUCTURE_COMPUTATION_CHECK_20260909.json)，根 `/structureComputationQualification` |
| A | [独立gate Stage A](evidence/SELECTIVE_CORRECTION_CHECK_20260909.json)，根 `/selectiveCorrectionQualification` |
| P(b,m) | A根 + `/basePointPanels/{b}/{m}`，将括号变量替换为下文实际键 |
| Q | [更早关联信用](evidence/ASSOCIATION_EVIDENCE_CREDIT_KILL_20260908.json)，仅为旧阶段参照 |

R/S都是 \(\alpha=k\exp\theta\) 的无历史系数门控式，不是当前共享gate。R的对称几何切片、S的历史关联、A的暖候选面板也不能混作同一任务配置。旧Q的可塑载体更早，进入它仅用于理解演变，不借其数值支撑当前式。原始快照内的旧执行文字及判定标签保持原状，新导读不把它们转化为行动指令。

## 2. 关键证据与支持范围

| 版本、口径 | 数字与精确位置 | 能支持 / 不能支持 |
|---|---|---|
| R；\(\epsilon=0\)，六个归一化梯度的有效秩 | R `/tonicGeometry/effectiveRank`=1；`/axialGeometry/effectiveRank`=1.6680217887 | 改可塑载体改变局部写入几何；该切片 \(s=0\)，不证明胞体贡献或历史保持。 |
| S；两历史基点 | S `/validationBlocks/0/paperPrediction/c`=.9648481205，`/r`=.9897537322；原向量在同对象 `/a`、`/b` | 正向净变化的局部条件成立，但两方向仍高度相关。 |
| S；主步长，两笔原始变化 | S `/metrics/axialFull/main/H0`：B=3.4041131025e−9，D=8.5657925197e−11，keepRatio=.02516306674；同路径H1：B=3.4302684027e−9，D=1.5449352908e−10 | 有证据方向净变化，但量很小；不是地点正确率。 |
| S；四笔交错、等参数写入范数 | S `/metrics/{mode}/rhoKeep/min`；mode依次为axialFull、axialDetached、tonic：.07027563238、.07022987050、.000170067915 | 在此预算/序列的比较；full减detached仅约4.5762e−5，不是稳健实用优势。 |
| A；冷候选 | A `/basePointSummary/cold_candidates/S_cross`=.01978055397 | 低绝对交叉可数值复现；不是所有候选起态成立。 |
| A；暖候选两种分配 | A `/basePointSummary/warm_candidates_H0_H1/S_cross`=.53979175799；换成 `warm_candidates_H1_H0` 得.14316644510 | 某暖条件互扰明显，不能只报冷基点。 |
| A；暖H0/H1、同正号、H0先写 | A `/units/8/eta/retention/H0`：HOverU=.53628843344，DOverU=.46371156656 | 自身修正存在但其他写入大量回撤；未达原定H/U≤.25、D/U≥.50。 |
| 当前共享gate；局部解析 | 见§5输入指针与复算式：.02114478599/.02418554150/.02387539402 | 仅是新式在 \(\theta=0\) 的条件导数预言，无新式有限写入结果。 |
| 真实RGB/HM3D | 无对应原始视觉结果文件 | 未验证，用户已暂缓；合成特征/命令不能作真实渲染证据。 |

S的“keepRatio”是 \(D/B\)，四笔表仍除以第一次自身B，不是累计自身写入。因此不能把约.0703叫作正确关联保持率，或与A的 \(D/U\)直接相比。A中U累加自身教学的证据方向变化，H累加其他写入的有害部分，D是净变化；交叉促进另记，不与H抵消。当前请求只是解释这些口径，不调整已冻结门槛。

## 3. 97.48%回撤、前向历史与胞体比例

从S的H0主更新原始对象：
\[
1-D/B=1-.02516306674244024=.9748369332575597.
\]
即第二笔回撤首笔修正约97.48%。不是“没有学习”，也不能用净变化为正掩盖强互扰。令S原始向量为 \(g_0,g_1\)，到达教学 \(\delta_0>0,\delta_1<0\)，
\[
c=\frac{g_0^Tg_1}{\|g_0\|\|g_1\|},\quad
r=\frac{\delta_0\|g_0\|}{|\delta_1|\|g_1\|}.
\]
两笔固定基点近似中 \(c<r<1/c\)使两段经历的净变化都顺教学方向；首笔线性回撤比例 \(c/r=.9748365569\)。此预言与上述实测接近，却并不意味着低交叉。根因不是符号或时序错误：S的共同成对参数方向主导，两支路贡献交叉内积约96.69%；可从 `paperPrediction/a,b` 的逐坐标乘积复算。历史已影响前向，但未转成足够不同的写入方向。

胞体贡献另核对S `/validationBlocks/0/fullDetached/fullGradient` 与 `/detachedGradient`：
\[
\frac{\|g_{\rm full}-g_{\rm detached}\|_2}{\|g_{\rm full}\|_2}
=.0031219338359.
\]
约0.312%是此L2复算值；同对象 `/fraction`=.003221397079使用源码的最大范数，约0.322%，不能标成L2。detached是同前向移除胞体误差响应的信用诊断，不是物理删除胞体。成对反对称读出使Schur分子容易抵消；\(\epsilon=0\)对称切片可有 \(\lambda_s=0\)，小不对称只弱打破抵消。解释见[合同§2](references/SELECTIVE_CORRECTION_GOAL_CONTRACT_20260909.md)，不能将新增均值池的价值转借给原误差通路。

## 4. 独立gate：暖候选的三路混合项

取 \(b=\mathrm{warm\_candidates\_H0\_H1}\)。原始梯度入口是P(b,m) `/rows/H0/gradient`、`/rows/H1/gradient`，其中m分别为 `full`、`commonOnlyOriginalR`、`differenceOnly`。在 \(\theta=0\) 令共同项 \(C_\nu\)、差分项 \(D_\nu\)，有精确分解 \(g_\nu=C_\nu+D_\nu\)。原向量内积给
\[
C_0^TC_1=4.8962606228\times10^{-6},\quad
D_0^TD_1=-1.2801377338\times10^{-5},
\]
\[
C_0^TD_1+D_0^TC_1=-6.8825058086\times10^{-7},\quad
g_0^Tg_1=-8.5933672964\times10^{-6}.
\]
最后一数也直接存于P(b,full) `/K/0/1`。并非共同项压过差分项，而是差分负交叉项及混合项改变了平衡；只改余弦符号没有解决选择性。

逐路入口是P(b,full) `/rows/{H0或H1}/routeCredits/{query或referenceA或referenceB}/localCredit`。以Q/A/B表示这三路向量，复算
\[
Q_0^TQ_1=-4.6326688e{-7},\quad Q_0^TA_1=-1.4070614e{-5},
\]
\[
B_0^TQ_1=-5.0544139e{-6},\quad A_0^TA_1=1.0330981e{-5}.
\]
完整九项都应求和，不能只挑最负项。query自身交叉很小，候选历史门控却把当前信用送往不同历史的坐标；共享参数下三路抵消失衡，是当前接口修订的因果依据，而不是模型已经证实修订有效。

同一单元A `/units/8/rho/retention/H0`给H/U=.37404990418、D/U=.62595009582，说明改等范数预算仍有互扰，不只是步长太大。eta原始U=1.5919645602e−9、H=8.5375218006e−10、D=7.3821238011e−10，均在 `/units/8/eta/retention/H0`；正自身写入没有消失。同负号及交换消费次序也在原24单元内，解读不能只选择反号互助单元。

## 5. 共享gate解析数怎样得到

对三基点b依次取cold_candidates、warm_candidates_H0_H1、warm_candidates_H1_H0，从P(b,commonOnlyOriginalR) `/rows/{H0或H1}/gradient`读 \(g_\nu^R\)，从P(b,full) `/rows/{H0或H1}/cells/query/w`读 \(w_\nu\)。在零参数前向相同的条件下：
\[
g_\nu^{\rm aligned}=w_\nu\odot g_\nu^R,\qquad
K_{\nu\omega}^{\rm aligned}=(g_\nu^{\rm aligned})^Tg_\omega^{\rm aligned}.
\]
按[理论篇](02_THEORY_AND_DERIVATION.md)定义计算S_cross，得到§2最后数值及[合同§8.4](references/SELECTIVE_CORRECTION_GOAL_CONTRACT_20260909.md)的矩阵。它不是把“旧独立gate总梯度”再乘query门控；使用的是原R总梯度，并仅在零参数同前向基点成立。对非零参数必须按新式重解三路。

## 6. 已覆盖与仍缺的检验

A `/correctnessChecks`数组有10项记录，但覆盖不完整。源码stageFiniteDifference约1822–1835行只对总分差作坐标差分，stageMargin默认容差为1e−12；P(b,full)各行 `/finiteDifferences`可逐项查误差。单独逐路FD及1e−14重算未覆盖，旧S的对应核对不能转借。源码event约855行系数域接受零值，而理论需要正有限指数；当前有限基点的互扰不据此归因于数值域问题。

旧Q/R/S回归分别保留118/132及S的检查与原始输出，详见[实现记录](references/SELECTIVE_CORRECTION_IMPLEMENTATION_20260909.md)。它们支持旧路径未被意外改变，不是当前共享gate证据，也不是能力测试总分。新式没有对应运行产物，不能把缺项写成已核验；下一步审阅应先规范数学与可判别证据，而不是自动续跑附件命令。
