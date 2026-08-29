# 交通流稳定性：可复现实验代码

本目录是《跟驰模型与宏观交通流模型的稳定性》的第一批开源配套代码。它不是另写的一套演示程序：书中的小算例、结果 JSON、插图和后续在线实验台均以这里的函数及实验清单为共同来源。

## 已实现内容

- `M01`：二阶局部欠阻尼响应，解析解与 RK4 逐点核对；
- `M02`：有限环道的离散波长效应；
- `M03`：原始 IDM 在稳定与失稳速度下的非线性环道积分；
- `M04`：不中断仿真状态，中途渐变开启前车加速度反馈；
- `M05`：用三个时间步长检查增长率收敛；
- `M06`：使用 IDM 平衡基本图的 LWR--Godunov 激波算例；
- `W01`：把三辆 IDM 环道的前两个 RK4 更新步完全展开；
- `W02`：把四网格 LWR 的界面通量和前两个有限体积更新步完全展开；
- `W03`：把三网格 ARZ 的守恒量、Rusanov 通量、松弛源项和前两个更新步完全展开；
- `D01`：NGSIM US-101 三分钟真实轨迹、连续车队和自然扰动增益分析；
- `D02`：NGSIM 八车队的单车／车队稳定性标定、留出递推验证和谱约束控制设计。

每个实验的入口和参考输出见 `experiments/manifest.json`。机器可读结果保存在 `results/`，书稿所用图保存在 `figures/`。

## 安装与运行

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python examples/run_mini_examples.py
python examples/worked_update_examples.py
python examples/ngsim_stability_application.py
python -m unittest discover -s tests -v
```

NGSIM 示例默认读取 `data/raw/us101_lane2_180s.csv`。若本地没有该文件，可以从美国交通部开放数据接口下载同一片段：

```bash
python examples/analyze_ngsim_us101.py --download
```

下载范围、字段、单位换算和筛选规则均写在脚本中。若只想复查书中的数值，不必保留 4 MB 原始切片，仓库已经包含两个小型派生表：八车队轨迹和 78 个跟驰片段的指标。

## 目录结构

```text
src/traffic_stability/     IDM、线性谱、环道 RK4 与 LWR 求解器
examples/                  M、W、D 系列可执行入口
experiments/manifest.json  书稿、代码与网站共用的实验清单
tests/                     解析极限和回归测试
data/derived/              可公开复查的小型派生数据
results/                   参考数值结果与逐步更新表
figures/                   由代码生成的图
```

## 复现约定

1. 所有速度、长度、时间和加速度计算使用 SI 单位；读取 NGSIM 后立即把英尺制字段转换为 SI。
2. 稳定性必须同时报告模型、平衡工况、边界条件、扰动、时间步长、拟合窗口和允许波数。
3. `M03`、`M05` 的增长率只在线性幅值窗口拟合；进入饱和后不再把包络斜率解释为线性特征值。
4. `D01` 的自然片段是观测证据，不等于对驾驶因果关系或全局线稳定性的证明。

## 数据来源与引用

NGSIM 轨迹由 U.S. Department of Transportation Federal Highway Administration 提供：

> U.S. Department of Transportation Federal Highway Administration. (2016). Next Generation Simulation (NGSIM) Vehicle Trajectories and Supporting Data. Dataset. DOI: 10.21949/1504477.

原始 NGSIM 文件不受本代码 MIT 许可证重新授权。使用者应遵守数据门户的引用说明和免责声明。`data/README.md` 给出本项目使用的精确查询。

## 许可证

代码采用 MIT License。书稿文字、排版和原始 NGSIM 数据不自动包含在该软件许可证中。正式公开前应把许可证中的项目署名替换为作者希望使用的姓名或机构名，并补充 `CITATION.cff`。
