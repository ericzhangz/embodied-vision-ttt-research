# 理论与推导：query条件的三路共享门控

> **历史方案，正文保留。** 当前候选已改为净耦合能量与关联诱发局部资格，请读 [最新简稿](01_CURRENT_IDEA_BRIEF.md)及[推导全集 §T](references/DERIVATION_PACKAGE.md)。本篇不是当前主公式，不将其预言转借给新方案。

## 1. 对象、状态与分类

总体研究对象是陌生静态室内导航中，后继证据修正跨视角关联并保持其他关联的在线过程。单次分差 \(f\)、混合证据损失 \(L_{\rm ev}\) 和局部Gram \(K\) 是分析对象，不替代地点正确性、连续保持或导航能力。当前公式的条件微分链可明确写出；选择性和结构额外价值仍是待证机制。

本篇用“定义/建模”“精确”“条件推导”“局部近似”“机制解释”区分层级。只整理已提出的算式，不新增能力定理。当前候选以[合同§8](references/SELECTIVE_CORRECTION_GOAL_CONTRACT_20260909.md)为准；实测归属见[证据链](03_EVIDENCE_CHAIN_AND_LIMITS.md)。

| 符号 | 含义及是否更新 |
|---|---|
| \(t;\xi\in\{t,A,B\};j=1,\ldots,4\) | 查询事件、三条读取路径、支路索引；A/B不是地点真值。 |
| \(x_\xi=(a,b),u_\xi,y_\xi\) | 冻结特征输入、实际执行动作、动作后已到达观测特征；目前输入是合成二维接口。 |
| \(Z_\xi^-;Z_\xi^+\) | 各自incoming与本token隐式末态；\(Z=(z_1,\ldots,z_4,s)\)。 |
| \(\theta\in\mathbb R^4\) | 在线可塑互易传递参数，直接参与描述子前向；不是输入突触权重。 |
| \(q,P,C_*,p,w\) | 历史活动、均值、冻结校准单位、固定接线极性、参数门控；\(w\)仅表示门控。 |
| \(h,k,\ell,\ell_s,\epsilon\) | token步长、基础传输、支路/胞体泄漏、固定输入调制强度。 |
| \(R,\hat D,f,\beta,J,M,\lambda\) | 固定读出、单位描述子、分差、状态偏导、流Jacobian、隐式响应矩阵、adjoint。 |
| \(\Lambda,\delta,g,\eta\) | 后继相对证据、标量教学、条件分差梯度、写入步长。 |

冻结视觉编码器、预处理、候选规则与读出；电位保存在线历史，只有 \(\theta\) 是此处TTT学习参数。没有另一个“慢参数创新”。学习/配对/选样不读取隐藏pose、depth、overlap、地点/实例标签、reward/success；动作相同只构成接口条件，不等于相同空间位置或物理位移。

## 2. 历史上下文进入同一前向

**定义/建模。** 当前固定接线
\[
V(x)=(a,-a,-b,b),\quad m(x)=(b,-b,a,-a),\quad p=(1,-1,-1,1).
\]
\(p\) 来自各输入的正负成对接线，不由历史编号或教学符号决定。只用查询incoming计算
\[
q_{t,j}=(z_{t,j}^--s_t^-)^2,\quad P_t=\tfrac14\sum_jq_{t,j},\quad
e_{t,j}=(q_{t,j}-P_t)/C_*,\quad w_{t,j}=1+p_je_{t,j}.
\]
校准时 \(\theta=0\)，最初两个实际处理的观测token之后得到 \(q^*,P^*\)，冻结
\[
C_*=\sqrt{\tfrac14\sum_j(q_j^*-P^*)^2}>0.
\]
这是前缀活动对比单位，不是每次查询标准差；不读取后继或通过调整尺度挑结果。零尺度不定义此门控，不发明补救常数。

三路保留自身 \(Z_\xi^-,x_\xi\)，仅共享系数上下文：
\[
\alpha_{t,j}=k e^{\theta_jw_{t,j}},\quad d_j^\xi=z_j^\xi-s^\xi,\quad
I_j^\xi=kd_j^\xi+\alpha_{t,j}(d_j^\xi)^3,
\]
\[
F_j^\xi=-L_j^\xi z_j^\xi-I_j^\xi+V_j(x_\xi),\quad
L_j^\xi=\ell+\epsilon m_j(x_\xi),\quad
F_s^\xi=-\ell_s s^\xi+\sum_jI_j^\xi,
\]
\[
Z_\xi^+=Z_\xi^-+hF_\theta(Z_\xi^+,x_\xi;w_t).
\tag{1}
\]
条件为 \(h,k,\ell_s,L_j^\xi>0\)，有限 \(\theta,w\) 和正有限 \(\alpha\)。负 \(w\) 允许参数作用反向，不代表负电导。计算中的指数溢出/下溢需要显式界定；泄漏和电流项耗散不排除输入供能。隐式token是离散模型选择，不是旧连续ODE或RK4轨迹的精确替代。

**分支边界。** 旧R/S取 \(w=1\)；已测独立gate取每路 \(w^\xi=w(Z_\xi^-)\)；当前式(1)取同一 \(w_t\)。三者在 \(\theta=0\) 的互易前向相同，但参数导数不相同。纯差分 \(w=pe\) 是已撤回的分析分支。不要混用这些式子的信用或有限写入结果。

## 3. 归一化读取、Jacobian与条件梯度

**定义。**
\[
R=\tfrac12\begin{pmatrix}1&-1&0&0&0\\0&0&-1&1&0\end{pmatrix},\quad
r_\xi=RZ_\xi^+,\quad n_\xi=\|r_\xi\|>0,\quad \hat D_\xi=r_\xi/n_\xi,
\]
\[
f=\hat D_t^T(\hat D_A-\hat D_B).
\]
三路用同一 \(\theta\) 版本及查询上下文，省略描述子的条件下标 \( |t \)。零读出未定义。

**精确读出微分。** 记 \(T_\xi=(I_2-\hat D_\xi\hat D_\xi^T)/n_\xi\)，则
\[
\beta^t=R^TT_t(\hat D_A-\hat D_B),\quad
\beta^A=R^TT_A\hat D_t,\quad
\beta^B=-R^TT_B\hat D_t.
\tag{2}
\]
B负号和候选归一化不能省略。以下对单路暂省 \(\xi\)，所有incoming、观测、\(w_t,C_*\) 在偏导中固定。

**条件推导。** 在隐式末态令
\[
G_j=k+3\alpha_jd_j^2,\quad c_j=hG_j,\quad
\chi_j=1+hL_j,\quad A_j=\chi_j+c_j.
\]
流的非零Jacobian为
\[
J_{jj}=-(L_j+G_j),\quad J_{js}=J_{sj}=G_j,\quad
J_{ss}=-(\ell_s+\sum_jG_j).
\]
所以
\[
M=I_5-hJ=
\begin{pmatrix}\mathrm{diag}(A_j)&-c\\-c^T&1+h\ell_s+\sum_jc_j\end{pmatrix}.
\]
对任意非零向量 \(v=(v_1,\ldots,v_4,v_s)\)，
\[
v^TMv=\sum_j\chi_jv_j^2+(1+h\ell_s)v_s^2+\sum_jc_j(v_j-v_s)^2>0.
\]
因此此条件响应可逆。它不是Newton有限迭代总会收敛的保证，也不是跨整个历史展开的Jacobian。

定义 \(M^T\lambda=\beta\)。消去支路得到星形Schur：
\[
\lambda_j=\frac{\beta_j+c_j\lambda_s}{A_j},\qquad
\lambda_s=\frac{\beta_s+\sum_j(c_j/A_j)\beta_j}
{1+h\ell_s+\sum_jc_j-\sum_jc_j^2/A_j}
=\frac{\beta_s+\sum_j(c_j/A_j)\beta_j}
{1+h\ell_s+\sum_jc_j\chi_j/A_j}.
\tag{3}
\]
这是一标量胞体求解加支路回代；每路仍需自身状态、Jacobian和右端，不能称所有误差都天然是同一个标量或求解免费。

从式(1)微分：
\[
M\partial_{\theta_j}Z^+=h\partial_{\theta_j}F,\quad
\partial_{\theta_j}F=\alpha_jw_{t,j}d_j^3(-\mathbf e_j+\mathbf e_s).
\]
结合(2)(3)，
\[
g_j^\xi=-h\alpha_{t,j}w_{t,j}(d_j^\xi)^3
(\lambda_j^\xi-\lambda_s^\xi),\quad
g_j=\sum_\xi g_j^\xi=\partial_{\theta_j}f,\quad
\Delta\theta_j=\eta\delta g_j.
\tag{4}
\]
\(\mathbf e_j\) 是五维单位向量，不是历史对比 \(e_{t,j}\)。localCredit就是 \(g_j^\xi\)。它包含原胞体误差响应，但该响应非零不等于有额外性能。

这里不沿旧 \(\theta\to Z^-\to w_t\) 反传，不对校准求导，不是BPTT；候选读取却仍对当前 \(\theta\) 完整求导。新参数下须重算状态、Jacobian和三路信用，不能只给旧梯度乘门控。对每路定义 \(b_j^\xi=g_j^\xi/w_{t,j}\) 会在零门控处不妥，故直接定义
\[
b_j^\xi=-h\alpha_{t,j}(d_j^\xi)^3(\lambda_j^\xi-\lambda_s^\xi),\qquad
g_j=w_{t,j}\sum_\xi b_j^\xi.
\tag{5}
\]

## 4. Gaussian相对教师与因果消费

**建模。** A/B在查询后继到达前已完成，执行动作与查询一致；连续接口原定锁定每命令最早两个合格完成记录，不以学习结果更换。当前算术仍是固定A/B，不是已经实现该选择器。令
\[
v_C=y_C-x_C,\quad r_t^{\rm obs}=y_t-x_t,\quad
\mu_C=x_t+v_C,\quad \mathcal L_C=\mathcal N(y_t;\mu_C,I_2).
\]
动作后继 \(y\) 来自冻结观测特征，不是上述可塑描述子。教师与三路 \(\theta\) 导数断开；若未来改为可塑教师，下面的梯度不再完整。

**精确恒等式。** 固定两似然，令 \(\pi=\sigma(f)\)、\(Q=\pi\mathcal L_A+(1-\pi)\mathcal L_B\)，则
\[
\Lambda=\tfrac12(\|r_t^{\rm obs}-v_B\|^2-\|r_t^{\rm obs}-v_A\|^2),\quad
\pi^+=\frac{\pi\mathcal L_A}{Q}=\sigma(f+\Lambda),
\]
\[
L_{\rm ev}=-\log Q,\quad
\partial_fL_{\rm ev}
=-\frac{\pi(1-\pi)(\mathcal L_A-\mathcal L_B)}Q
=\pi-\pi^+=-\delta.
\tag{6}
\]
所以式(4)是该单事件条件损失的梯度下降。相同候选增量给 \(\Lambda=\delta=0\)；绝对预测误差大也不必产生相对教学。合法观测来源不证明单位方差增量具有地点辨识性。

一次到达链：先选完成记录并冻结查询incoming、门控、版本、\(f,g,\mu_C\)；执行动作后 \(y_t\) 确已arrival，再评价旧预测、计算 \(\delta\)、核对版本并消费一次；写入后从既有live状态继续。多个延迟pending跨参数版本如何处理仍需明确，不能静默用陈旧资格冒充当前梯度。复评旧关联使用其旧query上下文，不用最新query偷换评价对象。

## 5. 局部选择性及归因边界

**局部近似。** 对两段经历 \(\nu,\omega\in\{0,1\}\) 的冻结条件分差，
\[
K_{\nu\omega}=g_\nu^Tg_\omega,\quad
\Delta f_\nu=\eta\delta_\omega K_{\nu\omega}+O(\|\Delta\theta\|^2).
\]
展开要求分差在邻域二阶光滑、描述子不接近零、步长足够小；余项未给统一上界。两笔序列第二笔须重算，基点 \(K\) 不是全程保持定理。
\[
c=K_{01}/\sqrt{K_{00}K_{11}},\quad
S_{\rm cross}=\max(|K_{01}|/K_{00},|K_{01}|/K_{11}).
\]
符号翻转不减少绝对互扰；共同及相反教师均须有可用方向。

**机制解释。** 成对共同轴 \((1,1,0,0)/\sqrt2,(0,0,1,1)/\sqrt2\) 与差分轴 \((1,-1,0,0)/\sqrt2,(0,0,-1,1)/\sqrt2\) 只是原四坐标的分解。共同信用与历史交互可能相消交叉内积并保留自身范数；偏离成对对称时有混合项，不能省去。单纯正交换基不改变Gram，作用来自前向内的上下文参数化。

系数层精确为 \(\log(\alpha/k)=\theta\odot w_t\)，属于普通乘性门控；这不等同于整个隐式动态可被任意FiLM或线性网络替换。query条件读取是实际工程语义变更：统一本次关联的参数上下文，不手工路由历史ID。\(P\)需额外上报/汇总/广播，候选需每query重编码；同门控简单对照要有同样信息和成本。结构的额外价值只能由匹配对照区分，不能从链式法则或神经命名直接推出。
