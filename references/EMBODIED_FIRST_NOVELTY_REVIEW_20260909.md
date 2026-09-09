# 具身视觉优先查新：后继证据驱动的历史条件区室写入

检索日期：2026-09-09（Asia/Shanghai）。用途：供主协调者更新当前简稿；这是一份机制对照，非新颖性认证。

## 1. 范围、来源与结论边界

使用非 ARIS 的 `research` 技能：完整读取其 SKILL.md，依用户约束由本代理直接调查，未派生代理、改模型、执行实验、clone/install 或写入 .aris。项目机制仅依据本轮只读的 [CURRENT_IDEA_BRIEF.md](D:/EV-TTT/CURRENT_IDEA_BRIEF.md) 和任务消息；未打开其所链接的推导、结果、私有稿件或最早附件。协调者提供的公开代码定位链接已重新读取，未把既有索引的摘要当作本轮证据。外部查询只使用公开论文名称和通用术语，未提交私有机制文本。

主体覆盖 2025–2026 年空间 TTT、在线导航适应、持续场景记忆、VPR 关联筛选和 SLAM 后到证据修正；保留 TF-VPR 这一较早但机制重要的近邻。十篇视觉核心工作与其作者代码构成主体，脑计算只补三篇。检索包括 arXiv、CVF、PMLR、ICLR/NeurIPS 论文集、ECCV/ICML 官方目录和作者 GitHub。使用的通用检索组包括 `spatial memory test-time`、`navigation test-time training`、`visual place recognition online self-supervised`、`loop closure delayed learning`，以及公开题名精确查询。没有声称穷尽截至当日所有论文。

证据记号：**P**＝已读原文相关方法/实验段落；**A**＝原始摘要/作者项目页；**C**＝已读固定 commit 的关键实现与调用处；**V**＝官方 venue 目录/论文集；**作者声明**＝尚未由官方目录独立确认。会议状态按实际取得的证据报告，不从 arXiv 月份或第三方列表推断录用。

本轮最直接的三项对照是 **NavMorph 的在线场景记忆改写、TF-VPR 的地点伪关联筛选、FFN 的下一帧到达后适应**。Spatial-TTT 是空间参数记忆的重要近邻，但其录播空间问答与本方案导航目标不同。本轮未找到这十篇中完整实现本方案全部机制的工作；这只支持“候选差异未证”，不支持“首次”或“独创”。

## 2. 本方案四项候选主张

固定任务：HM3D 静态陌生室内 navigation 的在线空间情境/跨视角地点关联。不是 manipulation、ObjectNav、对象中心表示、显式三维建模或空间 QA。学习/配对只允许 RGB、执行动作及确实可测自身状态；隐藏 pose/depth/overlap、地点标签和 reward/success 不作教师。

| 编号 | 可核查的候选主张 | 必须区分的事实 |
|---|---|---|
| C1 | 当前查询与历史 A/B 的地点关联，可以由已完成同动作记录和随后实际到达的观测提供有方向的修正信号 | 数据来源合法，不等于信号能辨别地点；预测动力学相似，不等于同地点 |
| C2 | 后继教学消费此前缓存的历史条件资格；同一当前输入在不同历史条件下可获得不同关联写入方向 | 固定 incoming 的条件偏导，不是全历史 BPTT，也不是已实现连续因果事件流 |
| C3 | 可塑互易区室电流同时承担历史状态演化与信用分配，局部电流因子和胞体响应组成三路条件梯度 | 结构参与精确导数，不等于结构具有独立性能价值；四支路和三次幂不是必要创新定义 |
| C4 | 这种结构化修正在交错写入中保持先前关联，并带来导航收益 | 目前只有合成算术；实用/稳健的关联保持、真实视觉和导航效果均未成立 |

机制对照所用公式保持原方案不变：

\[
I_j=k(z_j-s)+\alpha_j(z_j-s)^3,\quad \alpha_j=k\exp\theta_j,\qquad Z^+=Z^-+hF_\theta(Z^+,x).
\]

已完成同动作记录给出 \(v_C=y_C-x_C\)、\(\mu_C=x_t+v_C\)。实际后继到达后，

\[
\Lambda=\tfrac12\big(\|y-x-v_B\|^2-\|y-x-v_A\|^2\big),\quad
\delta=\sigma(f+\Lambda)-\sigma(f),\quad
\Delta\Theta=\eta\delta\nabla_\Theta f.
\]

其中 \(f\) 是同版本 query/A/B 的归一化分差。局部导数载体含 \(-h\alpha_j(z_j-s)^3(\lambda_j-\lambda_s)\)，资格延迟后只消费一次。以上来自本地简稿；本轮没有复验实现或数值结果。

## 3. 十篇具身/视觉核心近邻

### E1. NavMorph — 在线场景记忆改写的直接近邻

**题名/日期/状态：** *NavMorph: A Self-Evolving World Model for Vision-and-Language Navigation in Continuous Environments*；2025-06-30 首发；ICCV 2025，由 [CVF 正式论文](https://openaccess.thecvf.com/content/ICCV2025/papers/Yao_NavMorph_A_Self-Evolving_World_Model_for_Vision-and-Language_Navigation_in_Continuous_ICCV_2025_paper.pdf) 核实。证据 P/C/V；[原文 §3.1–3.3、补充 §9.2](https://arxiv.org/html/2506.23468v1)。

任务是 MP3D 环境中的 R2R-CE/RxR-CE 连续导航。CEM 在训练及在线测试中检索场景上下文，并以前向加权更新改写记忆项；预训练世界模型另用视觉特征重建、动作监督和先验/后验匹配。不能把预训练预测损失写成测试时的误差梯度教师。

**重叠：** 在线场景历史参与当前表示，测试中持续修改可复用记忆。**差异：** 已读 CEM 路径按当前场景相似性写值，没有等待实际后继来重判历史 A/B 关联，也不改写区室传输参数。其语言、位置表示和导航图接口与当前约束不同。**缺口：** 本方案尚未证明比这种简单记忆更新更准确地纠正地点混淆；代码也不足以支持 NavMorph 任意跨 episode 不重置的强说法。作者 README 的 Online Evaluation 同样提及 pseudo iterative demonstrator；本轮未追完整监督来源及位姿依赖，不能据 CEM 写入函数认证其评测全链 GT-free，或与我们的 RGB/action-only 权限等价（见 C-N）。

### E2. TF-VPR — 地点伪关联构造与筛选的直接近邻

**题名/日期/状态：** *Self-Supervised Place Recognition by Refining Temporal and Featural Pseudo Labels from Panoramic Data*；2022-08-19 首发、2024-11-20 v3；作者页标 RA-L 2024，arXiv 链接正式 DOI `10.1109/LRA.2024.3495584`。证据 A/C；[arXiv 记录](https://arxiv.org/abs/2208.09315)、[作者项目页](https://ai4ce.github.io/TF-VPR/)。这是 ai4ce 的 TF-VPR，绝非 2026 年同名 training-free benchmark。

以时间邻域启动伪正例，迭代表征学习、特征邻域扩张和候选验证；包含 Habitat-Sim RGB 与真实全景 VPR。实际范式是对已采集数据反复训练/挖掘，不能称为单遍导航闭环 TTT。

**重叠：** 不依靠地点标签的跨视角配对学习、当前表征驱动候选关联再筛选。**差异：** 本方案以实际后继的动作条件相对证据产生带符号的 A/B 信用；不是 epoch 级重挖。**关键限定：** 论文的“收缩”不应夸大成撤销既有 trusted pair；已读 RGB 实现保留旧 trusted 集合并添加通过验证的新项，见 §4。**缺口：** 本方案尚无真实跨视角伪教师校准证据。

### E3. Forget, Anticipate and Adapt / FFN — 后继观测决定历史适应的直接近邻

**题名/日期/状态：** *Forget, Anticipate and Adapt: Test Time Training for Long Videos*；2026-06-25 首发；[作者仓库](https://github.com/rajatmodi62/ffn) 与 [作者主页](https://rajatmodi62.github.io/) 标 ECCV 2026，本轮官方目录未取得同题独立确认，故保留“作者声明”。证据 P；[原文 §2.4–2.6、§3](https://arxiv.org/html/2606.26515v1)。

任务为长视频分割、分类和深度估计。下一帧实际到达后，与此前预测比较形成 surprise，决定是否适应；另缓存窗口退出帧的适应前特征，通过恢复损失主动遗忘。测试更新骨干/SSL 路径，下游任务头冻结。

**重叠：** 后到视觉证据触发对较早处理内容的参数适应，历史缓存参与更新。**差异：** 它控制适应时机与滑窗遗忘，没有同动作历史 A/B 地点身份的相对教师。**缺口：** 原文既称不看未来，又用 \(x_{t+1}\) 决定 \(x_t\) 的适应，合理解读需要至少一帧延迟；没有核实完整事件/输出调度，不能宣称零延迟。退出旧窗口的恢复也不是地点关联长期保持。

### E4. Spatial-TTT — 空间关联快权重的重要覆盖

**题名/日期/状态：** *Spatial-TTT: Streaming Visual-based Spatial Intelligence with Test-Time Training*；2026-03-12 首发；ECCV 2026，作者仓库与官方 [录用目录](https://eccv.ecva.net/Conferences/2026/AcceptedPapers) 同题记录核实，目录指向 poster/5426。证据 P/C/V；[原文 §3–4](https://arxiv.org/html/2603.12255v1)、[作者仓库](https://github.com/THU-SI/Spatial-TTT)。

SwiGLU 快权重吸收分块视觉 K/V，配合滑窗注意力和时空卷积；外层训练使用密集场景描述和空间 QA，测试内层写入来自视觉 token。评测为视频空间理解/流式空间感知；“route plan”问答不等于动作执行后的导航闭环。

**重叠：** 测试参数更新承载跨帧空间关联，不能再声称首次用 TTT 形成空间记忆。**差异：** 本方案以实际动作后继选择性修正历史候选，并学习区室传输；该工作的已读内层没有这种地点教师。**缺口：** “spatial-predictive”不应误读成实际下一动作观测误差；代码中的 V 是投影信号，详见 §4。

### E5. EmbodiedPlace — 跨视角候选修正，但不是在线参数学习

**题名/日期/状态：** *EmbodiedPlace: Learning Mixture-of-Features with Embodied Constraints for Visual Place Recognition*；2025-06-16 首发；本轮仅核实 arXiv，未取得正式 venue。证据 P；[原文 §3–5](https://arxiv.org/html/2506.13133v1)。

使用 GPS、时序、匹配或数据库自相似等约束选择候选邻居，以离线训练的 MoF 权重混合固定特征，再对查询候选重排；训练数据为 MSLS/Pitts-30k。并非陌生环境中持续执行梯度写入。

**重叠：** 通过相关历史视图修正地点候选表示。**差异：** 没有后继到达后的资格消费与参数信用，也不构成导航闭环。**缺口：** 不能把它全部说成依赖 GPS，原文明确有不同约束来源；也不能凭“具身”题名把它当成本方案相同协议的 TTT 基线。

### E6. Fast Spatial Memory / LaCET — 空间快权重保持的覆盖

**题名/日期/状态：** arXiv *Fast Spatial Memory with Elastic Test-Time Training*；2026-04-08 首发；ECCV 2026 的官方目录及 [作者页](https://fast-spatial-memory.github.io/) 使用 *Fast Spatial Memory with Scalable Elastic Test-Time Training*，作者一致。证据 P/V；[原文 §2–3](https://arxiv.org/html/2604.07350v1)。

面向 3D/4D 新视角/新时间渲染。每个 chunk 写快权重后，用在线重要性估计和锚点回拉约束漂移，锚点可采用 EMA；输入含相机/时间信息，外层用光度监督。既有显式 Gaussian 解码，也有直接渲染解码，不能一律说成显式场景表示。

**重叠：** 连续空间参数记忆及稳定性/可塑性问题。**差异：** 它按参数重要性稳定更新，不处理地点 A/B 后继信用。**缺口：** 本方案不能以极小正向净变化宣称已解决该类遗忘，更不能把重建指标与地点保持或导航收益直接相比。

### E7. FeedTTA — 导航闭环测试时学习已经存在

**题名/日期/状态：** *Test-Time Adaptation for Online Vision-Language Navigation with Feedback-based Reinforcement Learning*；PMLR 正式出版日期 2025-07-13 至 19，ICML 2025。证据 P/V；[正式条目及论文](https://proceedings.mlr.press/v267/kim25ad.html)。

导航 episode 完成后收到成功/失败二值反馈，用反馈 RL 更新策略，并正则梯度以平衡适应与稳定；评测包括 REVERIE。**重叠：** 陌生环境、动作闭环、延迟教学、测试适应。**差异：** success/failure 是本方案明确禁止进入学习的教师；它优化导航策略，不是地点关联传输。**缺口：** 本方案可以写“未使用成功反馈”，但无反馈不自动等于更好、更类脑或首个无监督导航适应。

### E8. IDEA — 近期跨环境适应经验积累

**题名/日期/状态：** *Turning Adaptation into Assets: Cross-Domain Bridging for Online Vision-Language Navigation*；2026-05-22 首发；ICML 2026 官方 [Downloads 目录](https://icml.cc/Downloads/2026) 取得同题 Poster 记录，指向 `/virtual/2026/poster/63184`；未读取最终论文集版。证据 P/V；[原文 §4](https://arxiv.org/html/2605.23257v1)。

冻结策略主体，通过 Fisher 加权的源/目标特征矩匹配优化 soft prompts，把 prompts、域坐标和质量指标保存为资产库。已有资产足够时组合使用，否则适应并写新资产；预算满后近邻合并。评测为 R2R/REVERIE/R2R-CE。

**重叠：** 历史条件适应、持续写入和跨环境复用。**差异：** 教师是源统计对齐与策略自身敏感性，不是动作后继对地点候选的裁决。**缺口：** “training-free”只描述可复用资产的分支，不能概括整个方法；本方案的历史保持尚无与此层面相当的证据。

### E9. ABot-Explorer / SG-Memo — 主动探索中的在线空间记忆

**题名/日期/状态：** *Explore Like Humans: Autonomous Exploration with Online SG-Memo Construction for Embodied Agents*；2026-04-21 首发；本轮核实 arXiv，未取得正式 venue。证据 P；[原文 §III–IV](https://arxiv.org/html/2604.19034v1)。

VLM 从 RGB 检出门、楼梯等通行节点，维护图并选择探索子目标；包括 HM3D/MP3D 与现实示例。测试记忆更新含节点关联、位置加权、合并及语义投票；VLM 权重学习在 SFT 阶段。虽然称 RGB-only，节点投影仍使用相机标定、外参位姿及地面假设。

**重叠：** 静态陌生室内主动探索、持续跨视角空间记忆。**差异：** 更新显式语义拓扑图，不是测试时学习隐式关联电流。**缺口：** 不能把“仅 RGB 感知”扩写为“无需位姿或几何假设”，也不能将它的闭环实证转借给本方案。

### E10. MASt3R-SLAM — 后到跨视角证据修正历史估计

**题名/日期/状态：** *MASt3R-SLAM: Real-Time Dense SLAM with 3D Reconstruction Priors*；2024-12-16 首发，2025-06-02 v2；CVPR 2025 由 [CVF 正式论文](https://openaccess.thecvf.com/content/CVPR2025/papers/Murai_MASt3R-SLAM_Real-Time_Dense_SLAM_with_3D_Reconstruction_Priors_CVPR_2025_paper.pdf) 核实。证据 P/V；[原文 §3.3–3.6](https://arxiv.org/html/2412.12392v2)。

单目 RGB 定位与稠密建图；新关键帧触发历史图像检索、匹配验证、回环边添加及全局优化，前端也持续融合 pointmap。后端优化关键帧 Sim(3) 位姿；教师是几何残差和模型匹配先验，不是在线改训 MASt3R 权重。

**重叠：** 后来的视图可以改变对先前空间状态的估计。**差异：** 修正对象是显式几何/位姿，不是地点表征的可塑关联参数。**缺口：** “后到证据修正历史”作为宽泛主张已被覆盖；SLAM 闭环是回环一致性，不等同于导航控制闭环。

### 补充筛选，不扩展核心表

- [ATENA / Active Test-time Vision-Language Navigation](https://arxiv.org/abs/2506.06630)，2025-06-07；[NeurIPS 2025 正式论文集](https://proceedings.neurips.cc/paper_files/paper/2025/hash/3f510f82323c1293ae3e343893dd77b1-Abstract-Conference.html) 核实。用选择性人类 episode 反馈及自评、混合熵校准；是 FeedTTA 附近的重要反例，不能把它统称成完全无外部反馈。
- [TTT3R](https://arxiv.org/html/2509.26645v1)，2025-09-30：将 CUT3R 状态更新解释为测试时关联学习，以置信度门控改进长序列重建。其“状态作为 fast weights”的算法对应值得注意，但不是地点后继教师；本轮未独立核实最终 venue，不按模型名字断言全网测试反传。
- [MTU3D](https://arxiv.org/abs/2507.04047)，2025-07-05：RGB-D query 记忆、目标 grounding 与探索一体化，跨 HM3D-OVON/GOAT/A-EQA；主要是训练轨迹上的学习和测试记忆状态更新，且目标中心/深度协议偏离本任务，未深读代码。
- [BudVLN](https://arxiv.org/abs/2602.06356)，2026-02-06，与 [Phi-Nav](https://arxiv.org/abs/2607.01754)，2026-07-02：虽含 retrospective/hindsight，原始摘要明确涉及 geodesic oracle 或 expert-guided on-policy 监督；不能当作仅依实际后继的无标签测试更新。
- [Scal3R](https://arxiv.org/abs/2604.08542) 与 [RetrieveVGGT](https://arxiv.org/abs/2605.09644) 为 2026 重建/检索状态方向；本轮仅筛查摘要层，未用于判定本方案精确机制被覆盖。AnyLoc 仅保留为静态无训练 VPR 参照，本轮未重复读其代码。

## 4. 本轮重新取得的固定源码证据

全部为公开 GitHub 网页/raw/API 的只读内容，没有执行作者代码。行号按 raw 文件实际换行复核；网页抽取行号可能不同。以下三项是代码支持的判断；其他论文不借用这个证据等级。

### C-N：NavMorph 的 eval 路径确实写记忆

固定 commit：`7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657`，通过作者仓库 commits/main API 取得。

- [ss_trainer_ETP.py L635](https://github.com/Feliciaxyao/NavMorph/blob/7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657/vlnce_baselines/ss_trainer_ETP.py#L635)：`_eval_checkpoint` 有 `torch.no_grad()`；仅此不能证明没有手写状态更新。
- [同文件 L1389](https://github.com/Feliciaxyao/NavMorph/blob/7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657/vlnce_baselines/ss_trainer_ETP.py#L1389)：`mode == 'eval'` 中将视觉和位置 embedding 拼接，先 `memory_vft_pos.push(...)`，再检索增强当前视觉表示。更新值来自当前观测表示，不是后继预测残差。
- [memory.py L153](https://github.com/Feliciaxyao/NavMorph/blob/7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657/utils_p/memory.py#L153)：`Memory_vft.push` 在满容量时按 key 余弦相似度取 5 个近邻，用 `alpha * old + (1-alpha) * new` 改 value；未满则插入。结构是有限关联记忆的前向改写。
- [加载处 L1002](https://github.com/Feliciaxyao/NavMorph/blob/7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657/vlnce_baselines/ss_trainer_ETP.py#L1002) 会读取预收集 CEM 文件。因此已证实测试 rollout 内更新，但未追完整 reset 生命周期，不宣称跨所有 episode 无限延续。
- [README Online Evaluation L109–113](https://github.com/Feliciaxyao/NavMorph/blob/7fb0f4aa44ff3e042335d8de0d1b2197c6f1f657/README.md#L109) 在评测说明中也写 pseudo iterative demonstrator，随后给出 `main.bash eval`。本轮只确认上述 eval 调用、写入对象及所读输入来源，未审清 demonstrator、监督与位姿依赖的完整传递链；不提供 GT-free 或 RGB/action-only 等价保证，也不反向断言此说明足以证明测试时使用 GT 教师。

### C-T：TF-VPR 的 epoch 挖掘、验证与旧关联保留

固定 commit：`83a4fbf59c1e51c0024c97943788db82fca9da7d`。

- [train_netvlad_RGB_ours.py](https://github.com/ai4ce/TF-VPR/blob/83a4fbf59c1e51c0024c97943788db82fca9da7d/train_netvlad_RGB_ours.py#L221)：已读 epoch 循环、`train_one_epoch`、数据库重编码、全库近邻及 `Compute_positive` 调用。学习对象是 ImageNetVlad 参数和候选集合；时间尺度为重复 epoch，不是等某个动作后继才消费一笔资格。
- [Verification_RGB.py L415](https://github.com/ai4ce/TF-VPR/blob/83a4fbf59c1e51c0024c97943788db82fca9da7d/models/Verification_RGB.py#L415)：`Verify_image` 调 SIFT、`matcher(..., threshold=0.5)`，匹配数至少 28 才保留。
- [matcher L496](https://github.com/ai4ce/TF-VPR/blob/83a4fbf59c1e51c0024c97943788db82fca9da7d/models/Verification_RGB.py#L496) 是 BF kNN 的描述子距离比率测试；已读函数没有 RANSAC/基础矩阵验证。文件名不能替代具体几何证据。
- [Compute_positive L565](https://github.com/ai4ce/TF-VPR/blob/83a4fbf59c1e51c0024c97943788db82fca9da7d/models/Verification_RGB.py#L565)：先排除旧 trusted，再验证新候选；没有新项则保留旧集，有则并入。由此不能把该 RGB 版本描述为“后到证据主动撤销已接受错误关联”。论文层面的候选收缩和实现中撤销旧标签须分开。

### C-S：Spatial-TTT 的实际写入信号与时序

固定 commit：`e2e33a62b6f92c33b7e24ff042be9737d05c5bdb`。

- [ttt_operation.py L71](https://github.com/THU-SI/Spatial-TTT/blob/e2e33a62b6f92c33b7e24ff042be9737d05c5bdb/qwen-vl-finetune/models/ttt_operation.py#L71)：`block_causal_lact_swiglu` 明确实现 Apply-then-Update，快权重为 `w0/w1/w2`。已读手写导数、可选 Muon/momentum 与范数恢复，是真实内层参数写入。
- [causal_swa_lact.py L921](https://github.com/THU-SI/Spatial-TTT/blob/e2e33a62b6f92c33b7e24ff042be9737d05c5bdb/qwen-vl-finetune/models/causal_swa_lact.py#L921)：decode 把 K/V 累积到 pending cache，满 chunk 才手写更新 `state.w0/w1/w2`；当前 query 的输出在该写入之前计算。
- [卷积定义 L382](https://github.com/THU-SI/Spatial-TTT/blob/e2e33a62b6f92c33b7e24ff042be9737d05c5bdb/qwen-vl-finetune/models/causal_swa_lact.py#L382) 与 [forward L1047](https://github.com/THU-SI/Spatial-TTT/blob/e2e33a62b6f92c33b7e24ff042be9737d05c5bdb/qwen-vl-finetune/models/causal_swa_lact.py#L1047)：深度可分离 Conv3d 聚合视觉 Q/K/V。已读更新直接使用 V 形成导数，没有显式 `prediction - observed_successor` 残差；因此本轮不称其为下一动作观测的 MSE 教师。
- 局部块的先读后写，不足以证明整条视频处理链具有本方案要求的逐动作因果性；3D 卷积窗口、视频 prefill 和块边界还需联合审查。本轮既不据此指控泄漏，也不将视频 QA 接口等同于连续机器人事件流。

FFN 额外读取了公开树和部分 trainer/model 文件，commit 为 `b53b99c97bfd9f164621086489ac4487c881dcf4`；尚未追到 MRM/AWA 与输出延迟的完整调用链，故 **FFN 的精确后继时序仍是论文证据，不计入上述三项完整代码证据**。

## 5. 三篇脑计算辅助对照

### B1. DendriCL：结构可以承载已有学习算法，但需可检验对应

*Dendritic In-Context Learning in a Single-Layer Spiking Neural Network*；2026-07-02，arXiv，未核实正式 venue。证据 P；[原文 §3.1–3.3](https://arxiv.org/html/2607.02283)、[日期/作者](https://arxiv.org/abs/2607.02283)。

顶树突状态实现 leaky online LMS，测试时突触参数冻结；训练用 BPTT，上下文提供带标签样本，任务为 Garg-style ICL。它说明“等价于已有优化”不自动否定结构贡献：算法状态落在何处、如何演化仍可检验。与本方案不同，DendriCL 改膜状态，本方案还改区室传输参数；二者不能因类脑词汇而互相替代。它不证明无标签地点教师有效，也不为本方案的胞体必要性背书。

### B2. CET：延迟资格与交错事件信用已经被明确研究

*Learning From the Past with Cascading Eligibility Traces*；2025-06-17；**ICLR 2026 正式发表**，由 [论文集](https://proceedings.iclr.cc/paper_files/paper/2026/hash/e647dad9086b5a4cc136e1d1926cc172-Abstract-Conference.html) 核实。证据 P/V；[原文 §3–4](https://arxiv.org/html/2506.14598v1)。

级联状态空间资格核把信用集中到延迟对应的历史事件，减少普通指数迹混合；实验含带标签视觉分类和有 reward 的 RL。覆盖“资格先保存、教学后到”的一般原则。当前显式缓存并一次性消费资格与 CET 的动态延迟核不同，但只是不同实现，还不能称更准确的长期信用；CET 也没有提供当前 A/B 地点教师。

### B3. Safaai 等：局部资格 × 胞体传来的区室误差已有精确分解

*Shunting Inhibition and Dendritic Branching Shape Local Credit Assignment*；2026-07-03 首发，2026-07-27 v2；arXiv，未核实正式 venue。证据 P；[本轮采用 v2 原文](https://arxiv.org/html/2607.03556v2)。

电导树网络把精确梯度分解为局部驱动/压差/电阻因子和沿树传输的胞体误差，区分精确传输与受限广播；任务为监督分类。覆盖“区室电导决定局部学习信用”的宽泛主张。其前馈树电导模型不等于本方案隐式互易星形动态、条件历史和延迟地点教师。原文也区分梯度重建与性能优势，不能由我们的梯度吻合直接推出胞体效用。另核实 [2607.24990](https://arxiv.org/abs/2607.24990) 存在，但未把其机制作为第四篇展开或借摘要作强结论。

## 6. 逐主张判断

| 主张 | 结论 | 根据与未闭合处 |
|---|---|---|
| C1：从测试视觉流自教空间/地点关联 | **已有覆盖；具体教师候选差异未证** | Spatial-TTT 覆盖参数关联；TF-VPR 覆盖无地点标签配对筛选；FFN 覆盖实际后继驱动适应。尚未核实有同动作 A/B 相对后继教师的完全同构近邻，但增量跨视角可迁移与地点判别正确性也未被本方案证明 |
| C2：历史条件资格延迟修正 | **已有覆盖；历史选择性候选差异未证** | CET 覆盖延迟信用；NavMorph 覆盖历史场景记忆更新。一次性缓存是因果记账约束，不是独立创新；从不同历史产生不同写入方向的实用选择性仍需证据 |
| C3：区室结构对应信用运算 | **候选差异未证** | Safaai 已给出电导局部资格与胞体传输分解，DendriCL 已明确结构内嵌在线算法。互易非线性流、三路条件响应与地点后继的具体对应尚可核查，不能仅按算式可化为梯度而否定，也不能按生物命名通过 |
| C4：交错保持与真实导航收益 | **候选差异未证** | FSM/IDEA 已研究持续记忆稳定性；当前只有合成算术，且强串扰。QA、重建、监督分类与 HM3D 地点学习是**范围不同**，不能把这种差异冒充已测优势 |

本轮没有看到完整同构先例，主要意味着这些近邻分布在不同更新对象和教师接口上；不能据此认定把部件组合起来就形成可发表结构创新。反过来，几何优化、梯度等价和已有资格迹，也不足以排除一种经验证的结构—运算对应。

## 7. 目前不能宣称什么；只保留一个差异继续核查

不能宣称：首个空间 TTT；首个无地点标签的跨视角关联学习；首个后到观测修正历史；首个局部电导/胞体信用机制；已解决保持/遗忘；稳健的胞体等预算优势；仅局部计算所以没有反向/全局响应成本；已完成 HM3D 或 RGB 导航实证。也不能声称 TF-VPR 所读版本会撤销旧 trusted pair，或 Spatial-TTT 的内层教师就是下一动作的真实视觉残差。

简稿报告的证据边界必须留在结论中：主更新证据方向净变化约 **8.57e-11、1.54e-10**；第二笔撤回首笔修正约 **97.48%**；full 对 detached 的等写入量差约 **4.58e-5**，主更新排序反转。胞体改动信用约 **0.312%** 说明参与，不说明效用。D/B 不是地点保持率。这些数值仅沿用简稿，不是本轮独立复现。

**唯一保留的候选差异：在同一持续空间关联算子中，实际动作后继提供的相对证据，能否通过历史条件化的互易区室响应，选择性改写此前混淆的地点关联。**

这里的结构—运算对应是：历史状态改变局部压差及微分电导；它们连同胞体误差响应决定该笔关联资格；实际后继只给该笔信用的标量方向/强度；写入再改变未来的同一传输算子。候选价值应落在这种对应是否产生可复用的历史选择性，而非四支路、三次幂、sigmoid 或 Schur 运算本身。

这句话仍须同时具备三类证据才能升级为结果：实际后继确实含有地点关联的判别信息；该信息在连续因果事件流中写回正确历史关系；区室响应对修正与保持有稳定且非微量的实用贡献。当前三项均未闭合。提出这些证据门槛不更换 benchmark、不增加候选方案，也不把尚未实现的连续事件驱动当作已有能力。

## 8. 给主协调者的七条结论

1. 最近三项对照应为 **NavMorph、TF-VPR、FFN**：分别对照在线场景记忆改写、地点伪配对、实际后继触发适应；Spatial-TTT 单列为参数空间记忆覆盖，避免主导全部叙事。
2. **NavMorph 的 eval 写入已从固定代码确认**；所读 CEM 路径属于近邻 value 加权更新，不能把预训练预测损失写成该路径的部署教师。README 的 Online Evaluation 也提及 pseudo iterative demonstrator；完整监督/位姿依赖未审清，不认证其与我们的数据权限等价。
3. **TF-VPR 比名称判断更接近地点学习，也比笼统“纠错”叙事更有限**：所读 RGB 代码筛新候选并保留旧 trusted，SIFT 比率/数量验证不等于 RANSAC，更不等于撤销已接受关系。
4. **FFN 阻止“首次用实际下一观测触发过去适应”的主张**；其时序至少需解释一帧延迟，完整代码调度本轮未确认。它的滑窗遗忘与地点长期保持范围不同。
5. **已有导航 TTA 和持续记忆不能遗漏**：FeedTTA 使用 episode 反馈；IDEA 存储并复用适应资产；ABot-Explorer 做在线图记忆并主动探索。不能概括成现有工作都只做录播 QA。
6. **脑计算只保留结构证据标准**：DendriCL 展示状态内嵌 LMS；CET 已正式发表于 ICLR 2026；Safaai 已覆盖局部资格/胞体传输分解。优化等价不自动否定结构创新，生物结构也不自动成立贡献。
7. 简稿只宜保留 §7 的**单一候选差异**，并同步保留微小效应、97.48% 回撤、胞体优势不稳和未有 RGB/HM3D 的限制；逐项结论用“已有覆盖 / 候选差异未证 / 范围不同”，不用主观分数或新颖性认证。

来源限制：ECCV/ICML 部分 poster 页面在浏览器抽取中返回 internal error，会议状态改由官方目录原始 HTML 和作者对应信息核实；FFN 仍只保留作者 venue 声明。无法打开某个页面或未追完调用链的情况，没有用第三方自动摘要替代原文结论。
