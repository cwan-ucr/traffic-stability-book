# NGSIM 数据说明

`D01` 使用美国交通部开放数据集 `8ect-6jqj` 中 US-101 南向车道 2 的 179.5 秒轨迹片段。

- 数据集 DOI：<https://doi.org/10.21949/1504477>
- API：<https://data.transportation.gov/resource/8ect-6jqj.csv>
- 地点筛选：`location="us-101"`
- 车道筛选：`lane_id="2"`
- 时间筛选：`global_time between 1118846979700 and 1118847159700`
- 访问日期：2026-08-28

原始字段中的 `local_y`、`v_length`、`v_vel`、`v_acc` 和 `space_headway` 分别按英尺、英尺、英尺每秒、英尺每二次方秒和英尺解释，进入分析后乘以 0.3048 转为 SI 单位。`global_time` 的相邻帧间隔为 100 ms。

## 派生文件

- `derived/us101_selected_platoon.csv`：在约第 100 秒处自动寻找最长前后车链，选取其中连续 20 秒的 8 辆车，按领车到上游跟车排序。
- `derived/us101_pair_episode_metrics.csv`：前车编号不变、采样连续且不少于 20 秒的跟驰片段；速度用 21 帧三阶 Savitzky--Golay 滤波后线性去趋势，再计算标准差增益及 0--5 秒范围内的最大相关滞后。

原始切片位于 `raw/`，默认被 `.gitignore` 排除；运行 `python examples/analyze_ngsim_us101.py --download` 可以从官方接口重新取得。派生文件只用于教学复核，不应替代官方原始数据。

