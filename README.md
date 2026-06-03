# 🛰 遥感数据预处理工作台

## 概述

这是一个集成了三个预处理模块的 Web 界面，用于遥感影像的批量处理：

```
大气校正 → 几何校正 → 干扰剔除
```

支持 **Landsat 8/9** 和 **Sentinel-2** 多光谱影像，支持 **.npy 格式**数据。

## 安装与运行

### 环境要求

- **Python 3.9+**（在 macOS / Python 3.9.6 上测试通过）
- 依赖 **rasterio 1.3.8**(读写 GeoTIFF/ENVI)。pip 安装的 rasterio wheel 已自带 GDAL，**无需单独安装 GDAL**；若使用源码编译版本，则需系统先装好 GDAL。
- 蚀变/InSAR 等功能依赖 NumPy、SciPy、Matplotlib、Pillow；可选的矿床类型 LLM 推理依赖 openai/anthropic 客户端(见下方环境变量)。

### 1. 获取代码

```bash
git clone git@github.com:KevinJH82/geo-analyser.git
cd geo-analyser
```

### 2. 创建虚拟环境并安装依赖(推荐)

```bash
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
```

> 不用虚拟环境也可以直接 `pip install -r requirements.txt`，但推荐隔离以避免污染系统环境。

### 3. 启动应用

```bash
python3 app.py
```

或运行启动脚本(等价)：

```bash
bash run.sh
```

服务默认监听 `0.0.0.0:5001`(`debug=True`)，因此**本机与同局域网设备**均可访问。

### 4. 打开浏览器

- 本机：**http://127.0.0.1:5001**
- 局域网其他设备：`http://<本机IP>:5001`

> ⚠️ `debug=True` 仅适合本地/内网开发使用，请勿直接暴露到公网。

### 5. 可选环境变量

以下变量均有默认值，按需在启动前 `export`(或写入 shell 配置)即可：

| 变量 | 用途 | 默认值 |
|------|------|--------|
| `DELIVERY_ROOT` | 交付项目数据根目录 | `/Volumes/大硬盘可劲用/DeepExplor/数据留存/交付数据` |
| `RESULTS_ROOT` | 蚀变分析结果落盘目录 | 应用目录下 `results/` |
| `MAX_RUNS_PER_DEPOSIT` | 每个矿床保留的历史结果数 | `30` |
| `DEEPSEEK_API_KEY` | 矿床类型 LLM 推理(未设置则自动跳过该功能) | 空 |
| `DEEPSEEK_MODEL` / `DEEPSEEK_BASE_URL` | LLM 模型名 / 接口地址 | `deepseek-chat` / `https://api.deepseek.com` |

示例：

```bash
export DELIVERY_ROOT="/path/to/交付数据"
export DEEPSEEK_API_KEY="sk-..."     # 仅在需要 LLM 辅助识别矿床类型时设置
python3 app.py
```

## 使用流程

### Step 1: 设置数据路径

- **输入目录**: 选择包含 `.npy` 文件的目录
  - 支持递归扫描子目录
  - 点击【扫描】按钮发现所有 `.npy` 文件
  
- **输出目录**: 设置处理结果的保存路径
  - 自动创建不存在的目录
  - 原始子目录结构会保留

### Step 2: 选择文件

在文件列表中勾选要处理的影像：

```
[✓] 2024/scene1/image.npy      (45.2 MB)
[✓] 2024/scene2/image.npy      (42.8 MB)
[ ] 2024/scene3/image.npy      (48.1 MB)
```

### Step 3: 配置处理参数

#### 📡 大气校正

| 参数 | 说明 | 默认值 |
|------|------|--------|
| 校正方法 | DOS 或 Simple | DOS |
| 太阳天顶角 | 太阳相对于天顶的角度（度） | 30.0° |
| 气溶胶光学厚度 | AOT550（仅 Simple 方法） | 0.1 |
| 地表高度 | 高程（km，仅 Simple 方法） | 0.0 |

**DOS 方法**: 基于暗目标消减，适合矿产勘探  
**Simple 方法**: 快速辐射校正，精度较低

#### 📍 几何校正

| 参数 | 说明 | 默认值 |
|------|------|--------|
| 校正方法 | Affine 或 Polynomial | Affine |
| 像元大小 | 输出像元分辨率（m） | 30.0 |
| 插值方式 | 最近邻 / 双线性 / 三次 | 双线性 |

**Affine**: 线性变换，适合小误差  
**Polynomial**: 高阶多项式，适合非线性畸变

#### 🧹 干扰剔除

| 参数 | 说明 | 默认值 |
|------|------|--------|
| 传感器 | Landsat 8/9 或 Sentinel-2 | Landsat 8/9 |
| NDVI 植被阈值 | 高于此值判定为植被 | 0.20 |
| 云蓝波段阈值 | 高于此值判定为云 | 0.25 |

### Step 4: 开始处理

点击【▶️ 开始处理】按钮，进度条和日志实时显示处理状态。

## 输入/输出格式

### 输入格式

- **文件格式**: `.npy`（NumPy 二进制格式）
- **数组格式**: `(bands, rows, cols)` 三维数组
- **数据类型**: 
  - 大气校正输入: `uint16` （DN 值，范围 0-65535）
  - 其他步骤: `float32` （反射率，范围 0-1）

### 输出格式

- **文件格式**: `.npy`
- **命名规则**: `{原文件名}_corrected.npy`
- **目录结构**: 完整保留输入目录结构

#### 示例

```
输入目录 (./input/)
├── 2024/
│   ├── scene1/
│   │   └── band_stack.npy
│   └── scene2/
│       └── band_stack.npy

输出目录 (./output/)
├── 2024/
│   ├── scene1/
│   │   └── band_stack_corrected.npy
│   └── scene2/
│       └── band_stack_corrected.npy
```

## 处理流程详解

### 1️⃣ 大气校正

**输入**: DN 值影像 (uint16)  
**输出**: 地表反射率 (float32, [0, 1])

移除大气干扰，将 DN 值转换为物理意义的反射率。

**DOS 方法**：
- 自动识别暗目标（水体、阴影）
- 通过暗像元 DN 值推断大气顶部辐射
- 逐波段校正

**Simple 方法**：
- 基于气溶胶光学厚度（AOT）的大气透过率模型
- 快速计算，精度相对较低

### 2️⃣ 几何校正

**输入**: 反射率影像 (float32)  
**输出**: 地理配准影像 (float32)

将影像配准到地理坐标系统。

**自动 GCP 生成**：
- 使用影像四角作为地面控制点
- 根据像元大小和投影参数计算地理坐标
- 多项式方法使用 9 个点覆盖整个影像

**Affine 变换**：简单、快速  
**Polynomial 变换**：处理复杂畸变（镜头畸变、轨道倾斜）

### 3️⃣ 干扰剔除

**输入**: 反射率影像 (float32, [0, 1])  
**输出**: 干净影像 (float32)，干扰区域为 NaN

识别并移除：
- 🌱 植被（NDVI > 阈值）
- 💧 水体（NDWI + MNDWI）
- 🏗️ 建筑（NDBI + NDVI）
- ☁️ 云（高蓝波段反射率）
- ❄️ 雪/冰（NDSI）

## 故障排除

### Q: "找不到 .npy 文件"

**A**: 
- 确认目录路径正确（支持 `~` 展开）
- 检查文件确实是 `.npy` 格式
- 确保有读取权限

### Q: "处理失败：形状不匹配"

**A**: 
- 确保输入影像是 `(bands, rows, cols)` 格式
- Landsat 8/9 通常是 7 波段
- Sentinel-2 通常是 11 波段

### Q: "scipy 未安装"

**A**: 
```bash
pip install scipy
```

### Q: "输出文件很大 / 包含 NaN"

**A**: 
- `.npy` 是无损压缩，干扰区域记录为 NaN
- 可以用 `np.nan_to_num()` 填充或用掩膜处理
- 在 GIS 软件中将 NaN 视为 nodata 值

## 技术架构

```
Frontend (HTML + CSS + JavaScript)
        ↓
Flask Web Server (app.py)
        ↓
Processing Pipeline
    ├── atmospheric_correction.py
    ├── geometric_correction.py
    └── interference_removal.py
        ↓
File I/O (.npy)
```

## 依赖

| 包 | 版本 | 用途 |
|----|------|------|
| Flask | 2.3+ | Web 框架 |
| NumPy | 1.24+ | 数组计算 |
| SciPy | 1.11+ | 插值和重采样 |
| Matplotlib | 3.7+ | 图像预览生成 |

## 输出日志示例

```
⏳ 读取文件: 2024/scene1/image.npy
【步骤 1】大气校正 (2024/scene1/image.npy)
【步骤 2】几何校正 (2024/scene1/image.npy)
【步骤 3】干扰剔除 (2024/scene1/image.npy)
保存结果: 2024/scene1/image.npy
✅ 完成: 2024/scene1/image.npy

📊 处理完成 - 成功: 1, 失败: 0
```

## 性能提示

- 单张 256×256 Landsat 8 影像处理时间：**~2-5 秒**
- 批量处理 10 张影像：**~20-50 秒**
- 内存占用：**~200-500 MB**（取决于影像大小）

处理时不要刷新浏览器，否则会中断。

## 常见用途

### 💎 矿产勘探

```
原始 Landsat 8 DN 值
  ↓ 大气校正 (DOS 方法推荐)
  ↓ 几何校正 (配准到地理坐标)
  ↓ 干扰剔除 (植被、水体移除)
  ↓ 矿物指数计算 (Fe、Al、Mg 等)
```

### 🌾 农业监测

```
原始 Sentinel-2 数据
  ↓ 大气校正
  ↓ 几何校正
  ↓ NDVI 计算
  ↓ 干旱评估、植被制图
```

## 许可

基于 NumPy、SciPy、Flask 等开源库构建。

---

**需要帮助？** 查看 geo-analyser 目录中的 README 或源代码注释。
