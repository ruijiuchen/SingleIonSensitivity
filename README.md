# SingleIonSensitivity — 逐帧投影与峰面积分析

## 概述

该程序读取一个 ROOT 文件中的二维直方图 (TH2F)，其 X 轴为频率 (MHz)、Y 轴为时间 (s)。
对于每个时间帧（每个 Y bin），将频谱沿 X 轴投影为一维直方图 (TH1F)，
在指定的频率范围内寻找峰，计算峰面积和本底面积，输出逐帧结果以及时间区间内的统计量。

## 输入数据

- **ROOT 文件**: `203Pt78_decay_BGSub/data/8243_PY82ch1_0242_trigger_11_2026-04-08T01-41-32.root`
- **TH2F 名称**: `h2d`（默认）
- **TH2F 结构**:
  - X 轴: 308.000000 ~ 311.999996 MHz, 209716 bins（~1.9e-5 MHz/bin ≈ 19 Hz/bin）
  - Y 轴: 0.000000 ~ 2.901253 s, 58 bins（~0.050 s/bin ≈ 50 ms/bin）
  - 数据: 背景减除后的计数（单位：任意强度）

## 计算流程

### 1. 时间窗口选择

用户通过 `--t1` 和 `--t2` 参数指定感兴趣的时间区间（单位：秒）。
程序将 `[t1, t2]` 映射到对应的 Y bin 范围 `[bin_t1, bin_t2]`。

### 2. 找峰（在平均谱上）

将所有选定时间帧的投影叠加为一张平均谱（`ProjectionX` 叠加），
在用户指定的频率范围 `[search_low, search_high]`（默认 309.732~309.736 MHz）内
扫描所有 bin，找到含量最高的 bin 的中心作为 **峰中心** `peak_center`。

### 3. 定义本底区域

本底区域中心取在峰中心右侧 `delta_f`（默认 2 kHz）处：
```
bg_center = peak_center + delta_f
```

### 4. 面积计算

对于每一帧的投影谱 `h1`，分别计算峰面积和本底面积。

#### 积分区间

- **峰积分区间**: `[peak_center - df/2, peak_center + df/2]`
- **本底积分区间**: `[bg_center - df/2, bg_center + df/2]`

其中 `df` 默认 0.0005 MHz（0.5 kHz）。

#### 面积值

```
area = sum_{b in [bin_low, bin_high]} h1.GetBinContent(b)
```

即积分区间内所有 bin 的含量直接加和。

#### 面积误差

以积分区间内各 bin 含量的**样本标准差**作为误差：
```
vals = [h1.GetBinContent(b) for b in range(bin_low, bin_high + 1)]
mean = sum(vals) / len(vals)
std = sqrt(sum((v - mean)^2 for v in vals) / (len(vals) - 1))
```

如果区间内只有 1 个 bin，误差为 0。

#### 净面积

```
net_area = peak_area - bg_area
net_error = sqrt(peak_err^2 + bg_err^2)
```

#### 相对误差

```
rel_err = err / area    (若 area = 0，则 rel_err = 0)
```

### 5. 时间区间平均

对所有选定帧的峰面积、本底面积、净面积分别计算：

| 统计量 | 公式 |
|--------|------|
| 均值 `mean` | `(1/n) * sum(v_i)` |
| 标准差 `std` | `sqrt(sum((v_i - mean)^2) / (n-1))` |
| 均值标准误 `mean_err` | `std / sqrt(n)` |

### 6. 输出

#### 终端打印

逐帧表（每帧一行）：时间中心、峰面积 ± 误差、本底面积 ± 误差、净面积 ± 误差、
峰积分频率边界（低边/高边）、本底积分频率边界（低边/高边）、
峰/本底/净的相对误差。

末尾统计行：均值、标准差、平均相对误差。

#### CSV 文件

`results_t{t1}_{t2}s_{文件名}.csv`，包含逐帧数据及末尾的 `# summary`、`# mean`、`# std` 行。

#### 帧图片

`frame_plots/frame_{编号}_t{时间}s_{文件名}.png`，每帧一张，包含：
- 蓝色直方图：该帧的 X 轴投影谱
- 绿色虚线（2条）：峰积分区间的左右边界
- 红色虚线（2条）：本底积分区间的左右边界
- 粉色点线：峰中心位置
- TPaveText 文字框：Peak/BG/Net 面积和误差

#### 面积演化图

`area_evolution_t{t1}_{t2}s_{文件名}.png`，包含三条 TGraphErrors 曲线：
- **Peak**（绿色空心圆 + 宽虚线）：峰面积随时间变化
- **BG**（红色实心方块）：本底面积随时间变化
- **Net**（蓝色实心菱形）：净面积随时间变化

#### ROOT 文件

`area_evolution_t{t1}_{t2}s_{文件名}.root`，保存了三条 TGraphErrors 和 TCanvas。

## 命令行参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--t1` | **必需** | 时间区间下限 (s) |
| `--t2` | **必需** | 时间区间上限 (s) |
| `--root-file` | data 目录下的指定文件 | ROOT 文件路径 |
| `--hist-name` | `h2d` | TH2F 名称 |
| `--search-low` | `309.732` | 找峰范围下限 (MHz) |
| `--search-high` | `309.736` | 找峰范围上限 (MHz) |
| `--df` | `0.0005` | 积分宽度 (MHz)，即 0.5 kHz |
| `--delta-f` | `0.002` | 峰中心到本底区域的偏移 (MHz)，即 2 kHz |
| `--plot-dir` | `frame_plots` | 帧图片输出子目录名 |

## 运行示例

```bash
conda activate base
python3 SingleIonSensitivity.py --t1 0 --t2 2.9
python3 SingleIonSensitivity.py --t1 0 --t2 1.35
python3 SingleIonSensitivity.py --t1 0 --t2 1.35 --df 0.0005 --delta-f 0.002
```

## 依赖

- ROOT (>= 6.32)，通过 conda 安装
- Python 3.8+
