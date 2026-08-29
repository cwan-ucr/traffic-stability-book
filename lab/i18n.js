(function () {
  "use strict";

  const pairs = [
    ["当前稳定", "Currently Stable"],
    ["当前失稳", "Currently Unstable"],
    ["当前已应用参数；精确解析判据", "Currently applied parameters; exact analytical criterion"],
    ["控制面板目标参数；精确解析判据", "Control-panel target parameters; exact analytical criterion"],
    ["实验台分区", "Laboratory sections"],
    ["目标稳定", "Target Stable"],
    ["目标失稳", "Target Unstable"],
    ["速度颜色：红=慢 · 蓝=快 · 平衡", "Speed color: red=slow · blue=fast · equilibrium"],
    ["仿真时间", "Simulation Time"],
    ["速度标准差", "Speed Standard Deviation"],
    ["当前前馈", "Current Feedforward"],
    ["平均速度", "Mean Speed"],
    ["净间距范围", "Net-Gap Range"],
    ["交通状态", "Traffic State"],
    ["重置", "Reset"],
    ["播放倍速只改变观看速度，不改变积分步长。环道放大模式会隐藏右侧面板，仿真仍继续。", "Playback speed changes only the viewing rate, not the integration step. Focus mode hides the side panel while the simulation continues."],
    ["启用前车加速度反馈", "Enable Leader-Acceleration Feedback"],
    ["应用在线干预", "Apply Online Intervention"],
    ["逐步关闭反馈", "Ramp Feedback Down"],
    ["当前理论", "Current Theory"],
    ["目标理论", "Target Theory"],
    ["当前 κ", "Current κ"],
    ["目标 κ", "target κ"],
    ["在 stop-and-go 已形成后再应用，可直接观察波动是否衰减。", "Apply it after stop-and-go waves form to observe whether the oscillation decays."],
    ["准备“不稳定 IDM”演示", "Prepare Unstable-IDM Demo"],
    ["600 s 时自动将 κ 平滑提高到", "Automatically ramp κ at 600 s to"],
    ["预设采用 50 veh/km，并给 0 号车一个 −2.5 m/s 的制动脉冲。手动流程：开始仿真 → 等待 σᵥ 明显增大 → 点击“应用在线干预”。自动流程会在 600 s 触发。", "The preset uses 50 veh/km and applies a −2.5 m/s braking pulse to vehicle 0. Manual workflow: start → wait for σᵥ to grow → apply online intervention. Automatic mode triggers at 600 s."],
    ["每次在线改变都会保留仿真时间、增益变化和过渡时长。", "Every online change records simulation time, gain change, and transition duration."],
    ["60 车开放车队，vₑ=8 m/s；头车施加 0.6 m/s、周期 48 s、持续 5 周期的正弦速度扰动。", "Open platoon of 60 vehicles at vₑ=8 m/s; the leader applies a 0.6 m/s sinusoidal disturbance with a 48 s period for five cycles."],
    ["HDV 末车振幅大于头车输入，说明扰动沿上游放大；协同 AV 把传播增益压到 1 以下并更快恢复。", "The last HDV oscillates more than the leader, indicating upstream amplification; cooperative AV control lowers the propagation gain below one and accelerates recovery."],
    ["轨迹颜色映射瞬时速度", "Trajectory color represents instantaneous speed"],
    ["末车扰动振幅", "Last-Vehicle Disturbance Amplitude"],
    ["恢复时刻", "Recovery Time"],
    ["Wu 模型扰动附加电耗", "Wu-Model Incremental Disturbance Energy"],
    ["协同 AV 不改变头车扰动，也几乎不改变平均速度；收益来自削弱逐车放大、提高 TTC，并把 Wu 模型估计的扰动附加电耗降低 89.4%。总电耗只下降 0.15%，因为大部分能量用于维持基准巡航。", "Cooperative AV control leaves the leader input and mean speed nearly unchanged. Its benefit comes from suppressing vehicle-to-vehicle amplification, increasing TTC, and reducing Wu-model incremental disturbance energy by 89.4%. Total energy falls only 0.15% because most energy sustains baseline cruising."],
    ["案例 B · 连续开放车流、低速瓶颈与固定空间 VSL", "Case B · Continuous Open Traffic, Slow-Vehicle Bottleneck, and Fixed VSL"],
    ["12 km 单车道高速公路，法定最高限速 120 km/h；取短车头时距 T=0.7 s、净间距 18 m，事件前均匀速度约 75 km/h。0–1200 s 持续到达（1087 辆新车进入）；慢车在 x≈6 km 处于 300–420 s 降至约 65 km/h。VSL 固定在 x=5–6 km，并严格在慢车事件于 t=300 s 发生后才开启。", "A 12 km single-lane motorway with a statutory 120 km/h limit uses T=0.7 s and an 18 m net gap, giving about 75 km/h before the event. Arrivals continue over 0–1200 s (1087 entering vehicles). The slow vehicle reduces to about 65 km/h at x≈6 km during 300–420 s. The VSL is fixed at x=5–6 km and turns on only after the event begins at t=300 s."],
    ["120 km/h 上限下，短时距工作点 Φ=−0.001415 s⁻²，处于长波失稳侧。t<300 s 时控制关闭；事件触发后才在 114–120 km/h 内执行由稳定裕度确定的 114 km/h 命令。该命令以 30 s 斜坡施加，再用连续流仿真事后评价排队、安全、延误和能耗。", "At the 120 km/h limit, the short-headway state has Φ=−0.001415 s⁻² and is long-wave unstable. Control is off for t<300 s. After the event, the 114 km/h command selected by the stability margin from 114–120 km/h is ramped in over 30 s. Queueing, safety, delay, and energy are evaluated afterward in simulation."],
    ["法定上限 120 km/h · 小 T 失稳", "Statutory Limit 120 km/h · Short-Headway Instability"],
    ["黑色虚线=向上游移动的排队尾部", "Black dashed line=upstream-moving queue tail"],
    ["最大稳定裕度 VSL 114 km/h", "Maximum-Margin VSL: 114 km/h"],
    ["绿色带=固定 VSL 区；黑色虚线=排队尾波", "Green band=fixed VSL zone; black dashed line=queue-tail wave"],
    ["在一公里控制区的运营约束 114–120 km/h 内最大化 Φ，得到 114 km/h；安全与运行指标不进入在线目标", "Maximizing Φ over the operational 114–120 km/h range in the one-kilometer control zone selects 114 km/h; safety and performance metrics are excluded from the online objective"],
    ["车辆 653 速度曲线", "Vehicle 653 Speed Time Series"],
    ["条件下按速度着色的道路位置时间轨迹图", " position-time trajectory diagram colored by speed"],
    ["指定车辆在无控制和稳定性控制下的速度时间曲线对比", "Selected-vehicle speed time-series comparison under no control and stability control"],
    ["Wh/车", "Wh/vehicle"],
    ["Wu 模型单位里程电耗", "Wu-Model Energy per Distance"],
    ["控制收益与代价", "Control Benefits and Tradeoffs"],
    ["法定上限仍是 120 km/h；真正导致失稳的是短时距 T=0.7 s 与高密度工作点。慢车在 t=300 s 触发扰动，VSL 此前保持关闭并从该时刻起响应。无控制排队尾部以 −2.69 m/s 向上游传播；事件触发型 114 km/h VSL 使尾波速度降到 −0.28 m/s、速度离差下降 53.9%、延误下降 30.1%。除 Φ 外均为事后读数。", "The statutory limit remains 120 km/h; instability is caused by the short T=0.7 s and high-density state. The slow vehicle triggers the disturbance at t=300 s, and the VSL remains off beforehand. Without control, the queue tail propagates upstream at −2.69 m/s. The event-triggered 114 km/h VSL reduces it to −0.28 m/s, lowers speed dispersion by 53.9%, and cuts delay by 30.1%. Every quantity except Φ is an ex-post simulation metric."],
    ["周长控制", "Ring-Length Control"],
    ["平均净间距", "mean net gap"],
    ["§4–5 基准与 IDM", "§4–5 Baseline and IDM"],
    ["扫描平衡速度，寻找扰动增长率变号的临界点", "Scan equilibrium speed to locate the sign change in disturbance growth rate"],
    ["测量：不稳定速度带", "Measurement: unstable speed band"],
    ["运行验证", "Run Validation"],
    ["运行后自动拟合线性增长段", "Automatically fit the linear-growth interval after the run"],
    ["仿真实测", "Measured in Simulation"],
    ["增长率由模态 m* 的复振幅 A(t) 拟合：Re λ = d ln|A|/dt，Im λ = d arg A/dt，波速 c₁ = −Im λ / k。 仅在扰动仍处线性阶段（|A| 未饱和）时有效。", "Growth rate is fitted from the complex amplitude A(t) of mode m*: Re λ=d ln|A|/dt, Im λ=d arg A/dt, and wave speed c₁=−Im λ/k. This is valid only while the disturbance remains linear and |A| is unsaturated."],
    ["min Ψ < 0 即失稳", "min Ψ < 0 indicates instability"],
    ["蓝=实测 红=理论斜率", "blue=measured; red=theoretical slope"],
    ["个模态", " modes"],
    ["全频段非负即线稳定", "Nonnegative over all frequencies implies string stability"],
    ["精确解析判据", "Exact Analytical Criterion"],
    ["长波门槛只取决于加速度核总量", "The long-wave threshold depends only on the total acceleration-kernel weight"],
    ["与延迟和分配距离无关", "and is independent of delay and allocation distance"],
    ["幅宽与波速缩放", "Amplitude-Width and Wave-Speed Scaling"],
    ["灰=初始剖面 · 红=所选时刻", "gray=initial profile · red=selected time"],
    ["基准参数", "Baseline Parameters"],
    ["同一 m=4 波形的波肩变尖，高次谐波能量相对基波由 0 增至 0.190。", "The shoulders of the same m=4 wave sharpen, while higher-harmonic energy relative to the fundamental rises from 0 to 0.190."],
    ["外形变陡说明非线性已介入，但不能单凭外形判成孤波。", "Steepening shows that nonlinear effects have entered, but shape alone does not establish a solitary wave."],
    ["稳定对照", "stable control"],
    ["0 车单次制动", "a single braking pulse on vehicle 0"],
    ["恢复定义为", "recovery is defined as"],
    ["连续 300 s", "for 300 consecutive seconds"],
    ["D01 · NGSIM US-101 真实轨迹", "D01 · NGSIM US-101 Real Trajectories"],
    ["辆车", " vehicles"],
    ["车道", "lanes"],
    ["每组连续时间不少于 20 s", "Each continuous episode lasts at least 20 s"],
    ["中位数；四分位区间", "Median; interquartile range"],
    ["中位数；bootstrap", "Median; bootstrap"],
    ["增益大于 1", "Gain above 1"],
    ["配对片段的中位数", "Median of paired episodes"],
    ["可复现处理链", "Reproducible Processing Pipeline"],
    ["对速度序列采用 21 帧三阶 Savitzky–Golay 平滑，并去除片段线性趋势；以跟驰车与前车速度标准差之比定义经验传递增益，时滞由互相关峰值估计。筛选仅保留连续不少于 20 s、前车速度标准差不少于 0.25 m/s 的同车道片段。", "Speed series are smoothed with a 21-frame cubic Savitzky–Golay filter and detrended within each episode. Empirical transfer gain is the follower-to-leader speed standard-deviation ratio; delay is estimated from the cross-correlation peak. Screening retains same-lane episodes lasting at least 20 s with leader speed standard deviation of at least 0.25 m/s."],
    ["样本中略多于一半的局部片段出现增益大于 1，但中位增益置信区间跨过 1。因此它支持“扰动可能逐车放大”的现象性证据，却不足以单独判定整段车流串行失稳；还需控制车道变换、共同外部激励与测量噪声。", "Slightly more than half of local episodes show gain above one, but the median-gain confidence interval crosses one. This supports the possibility of vehicle-to-vehicle amplification but cannot alone establish string instability for the entire traffic stream; lane changes, common forcing, and measurement noise must also be controlled."],
    ["官方数据 DOI", "Official Data DOI"],
    ["跟驰稳定性实验台", "Car-Following Stability Laboratory"],
    ["环道微观仿真 × 解析判据 · 理论与实测同屏对照", "Ring-road microsimulation × analytical criteria · theory and measurement side by side"],
    ["实时环道与在线控制", "Live Ring Road and Online Control"],
    ["实时环道", "Live Ring Road"],
    ["传播与探针", "Propagation and Probe"],
    ["理论与验证", "Theory and Validation"],
    ["非线性现象", "Nonlinear Phenomena"],
    ["应用案例", "Applications"],
    ["真实数据", "Empirical Data"],
    ["模型设置", "Model Settings"],
    ["目标稳定", "Target Stable"],
    ["目标失稳", "Target Unstable"],
    ["稳定", "Stable"],
    ["失稳", "Unstable"],
    ["环形道路实时仿真", "Live Ring-Road Simulation"],
    ["环形道路车辆运行状态，车辆颜色映射速度", "Vehicle states on the ring road; color represents speed"],
    ["速度颜色：红=慢 · 蓝=快", "Speed color: red=slow · blue=fast"],
    ["运行控制", "Run Control"],
    ["不中断仿真", "Simulation remains continuous"],
    ["暂停仿真", "Pause Simulation"],
    ["继续仿真", "Resume Simulation"],
    ["开始仿真", "Start Simulation"],
    ["退出专注", "Exit Focus"],
    ["放大环道", "Enlarge Ring Road"],
    ["播放倍速", "Playback Speed"],
    ["在线干预", "Online Intervention"],
    ["参数变化不会重置车辆状态", "Parameter changes do not reset vehicle states"],
    ["目标前馈增益", "Target Feedforward Gain"],
    ["平滑过渡时间", "Smooth Transition Time"],
    ["手动关闭", "Manual Off"],
    ["手动", "Manual"],
    ["自动", "Automatic"],
    ["当前参数已与目标一致", "Current parameters match the target"],
    ["目标参数尚未完全应用", "Target parameters are still being applied"],
    ["拥堵形成 → 控制消散演示", "Congestion Formation → Controlled Dissipation"],
    ["推荐入口", "Recommended Demo"],
    ["干预记录", "Intervention Log"],
    ["尚未干预", "No intervention yet"],
    ["扰动传播与车辆探针", "Disturbance Propagation and Vehicle Probe"],
    ["位置–时间轨迹图", "Position-Time Trajectories"],
    ["轨迹颜色映射速度", "Trajectory color represents speed"],
    ["车辆位置随时间变化的轨迹，轨迹颜色表示车辆速度", "Vehicle position-time trajectories colored by instantaneous speed"],
    ["速度时空图", "Speed Spatiotemporal Diagram"],
    ["间距时空图", "Gap Spatiotemporal Diagram"],
    ["各车辆速度随时间传播的时空图", "Propagation of vehicle speeds over time"],
    ["各车辆净间距随时间传播的时空图", "Propagation of net gaps over time"],
    ["速度与间距相对平衡值", "Speed and gap relative to equilibrium"],
    ["观测车辆", "Probe Vehicle"],
    ["指定车辆速度和净间距随时间变化曲线", "Speed and net-gap time series of the selected vehicle"],
    ["扰动包络", "Disturbance Envelope"],
    ["每辆车经历的最大速度偏差", "Maximum speed deviation experienced by each vehicle"],
    ["沿车辆编号的最大速度扰动包络", "Maximum speed-disturbance envelope along vehicle index"],
    ["当前读数", "Current Readings"],
    ["车辆编号", "Vehicle Index"],
    ["车号", "Vehicle Index"],
    ["车辆", "Vehicle"],
    ["时间", "Time"],
    ["速度", "Speed"],
    ["间距", "Gap"],
    ["道路位置", "Road Position"],
    ["平衡速度", "Equilibrium Speed"],
    ["长波裕度", "Long-Wave Margin"],
    ["长波稳定裕度", "Long-Wave Stability Margin"],
    ["AV 渗透率", "AV Penetration"],
    ["五类典型非线性稳定性现象", "Five Canonical Nonlinear Stability Phenomena"],
    ["选择非线性现象", "Select a nonlinear phenomenon"],
    ["波形陡化", "Wave Steepening"],
    ["非线性饱和", "Nonlinear Saturation"],
    ["孤立波与 kink", "Solitary Waves and Kinks"],
    ["有限幅触发与迟滞", "Finite-Amplitude Triggering and Hysteresis"],
    ["振幅依赖恢复", "Amplitude-Dependent Recovery"],
    ["高次谐波增长", "Higher-Harmonic Growth"],
    ["增长率随振幅下降", "Growth rate decreases with amplitude"],
    ["同一稳定工况、不同初幅", "Same stable condition, different initial amplitudes"],
    ["恢复时间", "Recovery Time"],
    ["初始剖面", "Initial Profile"],
    ["剖面时刻", "Profile Time"],
    ["高次谐波比", "Higher-Harmonic Ratio"],
    ["注入模态", "Injected Mode"],
    ["初始", "Initial"],
    ["末 300 s 均值", "Mean over final 300 s"],
    ["稳定对照末值", "Final stable-control value"],
    ["模板振幅", "Template Amplitude"],
    ["相对幅宽", "Relative Width"],
    ["尚未证实", "Not yet demonstrated"],
    ["完整 IDM 证据", "Full IDM Evidence"],
    ["约化方程模板", "Reduced-Equation Template"],
    ["触发证据；迟滞未证实", "Triggering observed; hysteresis not demonstrated"],
    ["交互控制", "Interactive Control"],
    ["不改变实时环道状态", "Does not modify the live ring-road state"],
    ["如何判读", "How to Read the Result"],
    ["现象 → 证据 → 限制", "Phenomenon → Evidence → Limitation"],
    ["统一实验口径", "Unified Experimental Protocol"],
    ["便于论文复现", "Designed for reproducibility"],
    ["稳定性分析的交通管理应用案例", "Traffic-Management Applications of Stability Analysis"],
    ["选择应用案例", "Select an application case"],
    ["案例 A · 开放车队协同稳定", "Case A · Cooperative Stabilization of an Open Platoon"],
    ["头车有限周期正弦扰动 · HDV 与多前车/后车协同 AV 对照", "Finite-cycle sinusoidal leader disturbance · HDV versus cooperative AV"],
    ["案例 B · 低速车辆与固定 VSL", "Case B · Slow Vehicle and Fixed VSL"],
    ["连续到达车流 · 事件点上游固定 1 km 限速区", "Continuous arrivals · fixed 1-km VSL zone upstream of the event"],
    ["案例 A · 开放车队的正弦扰动与协同 AV", "Case A · Sinusoidal Disturbance and Cooperative AVs in an Open Platoon"],
    ["多前车前馈 + 后车速度反馈", "Multi-predecessor feedforward + rear-vehicle speed feedback"],
    ["案例 B · 连续开放车流、低速瓶颈与固定空间 VSL", "Case B · Continuous Open Traffic, Slow-Vehicle Bottleneck, and Fixed VSL"],
    ["法定上限 120 km/h · 小 T 失稳", "Statutory limit 120 km/h · short-headway instability"],
    ["最大稳定裕度 VSL 114 km/h", "Maximum-margin VSL: 114 km/h"],
    ["应用数据未加载", "Application data not loaded"],
    ["完整 IDM · 开放边界 · 1200 s", "Full IDM · open boundary · 1200 s"],
    ["位置–时间轨迹 · 颜色=速度", "Position-time trajectories · color=speed"],
    ["每条线=一辆车的道路轨迹", "Each line is one vehicle trajectory"],
    ["黑色虚线=向上游移动的排队尾部", "Black dashed line=upstream-moving queue tail"],
    ["绿色带=固定 VSL 区；黑色虚线=排队尾波", "Green band=fixed VSL zone; black dashed line=queue-tail wave"],
    ["事件时刻以竖线标记", "Event times are marked by vertical lines"],
    ["与左图共用位置范围和速度色标", "Uses the same position range and speed scale as the left panel"],
    ["无控制", "No Control"],
    ["控制", "Control"],
    ["扰动开始", "Disturbance Starts"],
    ["扰动结束", "Disturbance Ends"],
    ["VSL 开始", "VSL Starts"],
    ["排队尾部", "Queue Tail"],
    ["形成段", "formation stage"],
    ["候选 VSL", "Candidate VSL"],
    ["推荐", "Selected"],
    ["120 km/h 上限下为何取 114 km/h", "Why 114 km/h Is Selected under a 120 km/h Limit"],
    ["允许范围内直接最大化稳定裕度", "Directly maximize the stability margin over the admissible range"],
    ["候选可变限速的长波稳定裕度与纯稳定性目标扫描", "Long-wave stability-margin scan over candidate VSL values"],
    ["红线 Φ=0 为中性边界；118–120 km/h 仍在失稳侧", "The red line Φ=0 is the neutral boundary; 118–120 km/h remains unstable"],
    ["指定车辆", "Selected Vehicle"],
    ["速度曲线", "Speed Time Series"],
    ["红=无控制 · 绿=稳定性控制", "Red=no control · green=stability control"],
    ["安全—效率—能耗读数", "Safety, Efficiency, and Energy Metrics"],
    ["事后仿真评价 · 不参与在线选速", "Ex-post simulation evaluation · excluded from online speed selection"],
    ["形成阶段排队尾波速度", "Formation-Stage Queue-Tail Wave Speed"],
    ["逆向传播减缓", "Upstream propagation reduced"],
    ["低速车辆峰值", "Peak Number of Low-Speed Vehicles"],
    ["峰值速度离差", "Peak Speed Standard Deviation"],
    ["最小 TTC", "Minimum TTC"],
    ["累计延误", "Cumulative Delay"],
    ["单位里程电耗", "Energy per Distance"],
    ["下降", "reduced by"],
    ["提高", "increased by"],
    ["提前", "earlier by"],
    ["管理含义", "Management Interpretation"],
    ["稳定性收益", "Stability Benefit"],
    ["控制收益与代价", "Control Benefits and Tradeoffs"],
    ["慢车在 t=300 s 触发扰动，VSL 此前保持关闭并从该时刻起响应。", "The slow vehicle triggers the disturbance at t=300 s; the VSL remains off beforehand and responds only from that instant."],
    ["事件触发后才在 114–120 km/h 内执行由稳定裕度确定的 114 km/h 命令。", "Only after the event does the controller apply the 114 km/h command selected from the 114–120 km/h range by the stability margin."],
    ["NGSIM 真实轨迹稳定性证据", "Stability Evidence from NGSIM Trajectories"],
    ["观测证据 · 不是对整体串行失稳的单独证明", "Observational evidence · not a standalone proof of network-wide string instability"],
    ["数据规模", "Data Scale"],
    ["有效跟驰片段", "Valid Car-Following Episodes"],
    ["时间车头时距", "Time Headway"],
    ["速度扰动增益", "Speed-Disturbance Gain"],
    ["响应时滞 / 相关", "Response Lag / Correlation"],
    ["估计方法", "Estimation Method"],
    ["怎样解释结果", "How to Interpret the Result"],
    ["统计不确定性必须保留", "Statistical uncertainty must be retained"],
    ["场景", "Scenarios"],
    ["基准 IDM", "Baseline IDM"],
    ["前车加速度前馈", "Leader-Acceleration Feedforward"],
    ["后车加速度", "Follower-Acceleration Feedback"],
    ["多前车 / 多后车", "Multiple Predecessors / Followers"],
    ["后向间距反馈", "Rear-Gap Feedback"],
    ["时间延迟", "Time Delay"],
    ["异质车流", "Heterogeneous Traffic"],
    ["混合交通 / AV 编组", "Mixed Traffic / AV Platooning"],
    ["环道", "Ring Road"],
    ["车辆数", "Number of Vehicles"],
    ["密度", "Density"],
    ["环道周长", "Ring Length"],
    ["车长", "Vehicle Length"],
    ["IDM 参数", "IDM Parameters"],
    ["期望速度", "Desired Speed"],
    ["安全时距", "Desired Time Headway"],
    ["最小间距", "Minimum Gap"],
    ["最大加速度", "Maximum Acceleration"],
    ["舒适减速度", "Comfortable Deceleration"],
    ["加速度指数", "Acceleration Exponent"],
    ["反馈与延迟", "Feedback and Delay"],
    ["前馈增益", "Feedforward Gain"],
    ["前向衰减", "Forward Decay"],
    ["后向衰减", "Backward Decay"],
    ["前视车辆数", "Number of Predecessors"],
    ["后视车辆数", "Number of Followers"],
    ["后向间距权重", "Rear-Gap Weight"],
    ["反应延迟", "Reaction Delay"],
    ["通信延迟", "Communication Delay"],
    ["异质与混合交通", "Heterogeneous and Mixed Traffic"],
    ["重型 / 保守车辆占比", "Heavy / Conservative Vehicle Share"],
    ["AV–HV 降级增益", "AV–HV Degraded Gain"],
    ["AV 编组数", "Number of AV Groups"],
    ["车辆排列方式", "Vehicle Arrangement"],
    ["成组", "Grouped"],
    ["随机", "Random"],
    ["连续", "Contiguous"],
    ["均匀", "Uniform"],
    ["扰动与数值", "Disturbance and Numerics"],
    ["初始扰动幅值", "Initial Disturbance Amplitude"],
    ["步长", "Time Step"],
    ["暂停", "Pause"],
    ["运行", "Run"],
    ["章节验证实验", "Chapter Validation Experiments"],
    ["选择验证实验", "Select a validation experiment"],
    ["中性稳定线", "Neutral Stability Boundary"],
    ["线性增长率", "Linear Growth Rate"],
    ["最不稳定波数", "Most Unstable Wavenumber"],
    ["逐车传递比", "Vehicle-to-Vehicle Transfer Ratio"],
    ["有限 N 效应", "Finite-N Effect"],
    ["异质排序效应", "Heterogeneous Ordering Effect"],
    ["AV 临界渗透率", "Critical AV Penetration"],
    ["延迟势阱", "Delay-Induced Stability Well"],
    ["后视闭合点", "Rear-Looking Closure Point"],
    ["走停波波速", "Stop-and-Go Wave Speed"],
    ["AV 渗透率–速度稳定性二维平面", "AV Penetration–Speed Stability Plane"],
    ["绿色稳定，红色失稳", "green=stable, red=unstable"],
    ["指向位置", "Pointer"],
    ["当前工作点", "Current Operating Point"],
    ["验证读数", "Validation Readings"],
    ["理论 vs 实测", "Theory vs Measurement"],
    ["拟合窗口", "Fitting Window"],
    ["解析预测", "Analytical Prediction"],
    ["相对误差", "Relative Error"],
    ["对应实验", "Experiment"],
    ["振荡频率", "Oscillation Frequency"],
    ["波速", "Wave Speed"],
    ["最不稳定模态", "Most Unstable Mode"],
    ["全速域不稳定带", "Unstable Speed Band"],
    ["模态振幅增长", "Modal-Amplitude Growth"],
    ["环道增长率谱", "Ring-Road Growth-Rate Spectrum"],
    ["混合环周传递裕度", "Mixed-Traffic Ring Transfer Margin"],
    ["频域判据", "Frequency-Domain Criterion"],
    ["目标工作点的解析量", "Analytical Quantities at the Target Operating Point"],
    ["实际", "Measured"],
    ["理论", "Theory"],
    ["文档基准", "Book Reference"],
    ["扫描结果", "Scan Result"],
    ["全速域稳定", "Stable over the full speed range"],
    ["按「运行」开始采样", "Press Run to start sampling"],
    ["请先运行约 30–60 s", "Run the simulation for about 30–60 s first"],
    ["等待线性拟合窗口", "Waiting for the linear fitting window"],
    ["待复核", "Requires Verification"],
    ["stop-and-go 已形成", "Stop-and-go waves formed"],
    ["扰动正在放大", "Disturbance is growing"],
    ["弱扰动", "Weak disturbance"],
    ["近似均匀流", "Approximately uniform flow"],
    ["当前", "Current"],
    ["目标", "Target"],
    ["辆", " vehicles"],
    ["号", ""],
    ["组", " groups"],
  ];

  const ordered = pairs.slice().sort((a, b) => b[0].length - a[0].length);
  const params = new URLSearchParams(window.location.search);
  const language = params.get("lang") || window.localStorage.getItem("tsl-language") || "zh";

  function translate(value) {
    let output = String(value == null ? "" : value);
    if (language !== "en" || !/[\u3400-\u9fff]/.test(output)) return output;
    for (const [source, target] of ordered) output = output.split(source).join(target);
    output = output
      .replace(/Current(Stable|Unstable)/g, "Current $1")
      .replace(/Target(Stable|Unstable)/g, "Target $1")
      .replace(/，/g, ", ")
      .replace(/。/g, ". ")
      .replace(/；/g, "; ")
      .replace(/：/g, ": ")
      .replace(/\s+([,.;:])/g, "$1")
      .replace(/\s{2,}/g, " ")
      .trim();
    return output;
  }

  function applyTranslations(root) {
    if (language !== "en") return;
    const scope = root && root.nodeType ? root : document.documentElement;
    if (scope.nodeType === Node.TEXT_NODE) {
      if (/[\u3400-\u9fff]/.test(scope.nodeValue || "")) scope.nodeValue = translate(scope.nodeValue);
      return;
    }
    if (scope.nodeType !== Node.ELEMENT_NODE && scope.nodeType !== Node.DOCUMENT_NODE) return;
    const walker = document.createTreeWalker(scope, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach((node) => {
      if (/[\u3400-\u9fff]/.test(node.nodeValue || "")) node.nodeValue = translate(node.nodeValue);
    });
    const elements = scope.querySelectorAll ? scope.querySelectorAll("[aria-label],[title],[placeholder]") : [];
    elements.forEach((element) => {
      ["aria-label", "title", "placeholder"].forEach((name) => {
        const value = element.getAttribute(name);
        if (value && /[\u3400-\u9fff]/.test(value)) element.setAttribute(name, translate(value));
      });
    });
  }

  window.TSL_I18N = {
    language,
    translate,
    setLanguage(next) {
      window.localStorage.setItem("tsl-language", next);
      const url = new URL(window.location.href);
      url.searchParams.set("lang", next);
      window.location.href = url.toString();
    },
  };

  document.documentElement.lang = language === "en" ? "en" : "zh-CN";

  if (language === "en" && window.CanvasRenderingContext2D) {
    const nativeFillText = window.CanvasRenderingContext2D.prototype.fillText;
    window.CanvasRenderingContext2D.prototype.fillText = function (text, x, y, maxWidth) {
      const translated = translate(text);
      if (arguments.length >= 4) return nativeFillText.call(this, translated, x, y, maxWidth);
      return nativeFillText.call(this, translated, x, y);
    };
  }

  const observer = new MutationObserver((records) => {
    if (language !== "en") return;
    records.forEach((record) => record.addedNodes.forEach((node) => applyTranslations(node)));
  });

  document.addEventListener("DOMContentLoaded", () => {
    applyTranslations(document.documentElement);
    observer.observe(document.body, { childList: true, subtree: true });
    if (language === "en") window.setInterval(() => applyTranslations(document.body), 500);
    document.querySelectorAll("[data-lang]").forEach((button) => {
      const active = button.getAttribute("data-lang") === language;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-pressed", active ? "true" : "false");
    });
  });
})();
