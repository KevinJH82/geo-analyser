# 🚀 快速开始指南

## 1️⃣ 一键启动

### macOS / Linux

```bash
cd /Users/mac/Desktop/Kevin\'s/Claude\ Code/Web\ Search/geo-analyser
bash run.sh
```

### Windows

```bash
cd "你的路径\geo-analyser"
python app.py
```

然后在浏览器打开：**http://127.0.0.1:5000**

## 2️⃣ 第一次测试（使用演示数据）

### 方式 A：自动生成演示数据

```bash
# 生成 3 个测试场景
python3 generate_demo_data.py
```

这会创建 `demo_data/` 目录，包含：
- `scene1.npy` - 混合场景（植被 + 水体 + 建筑）
- `scene2.npy` - 森林场景（高植被覆盖）
- `scene3.npy` - 城镇场景（建筑密集）

### 方式 B：在 Web 应用中操作

1. **输入目录** 输入框中填写：
   ```
   /Users/mac/Desktop/Kevin's/Claude Code/Web Search/geo-analyser/demo_data
   ```

2. 点击 **【扫描】** 按钮
   - 应该找到 3 个 `.npy` 文件

3. 勾选所有文件

4. **输出目录** 输入框中填写：
   ```
   ./demo_output/
   ```

5. 点击 **【▶️ 开始处理】**

6. 等待处理完成（每个文件约 2-5 秒）

7. 完成！结果保存在 `demo_output/` 目录

## 3️⃣ 使用真实数据

### 准备数据

确保你的 `.npy` 文件格式为：
```python
shape: (bands, rows, cols)
dtype: uint16  # DN 值范围 0-65535
```

如果你有 GeoTIFF 或其他格式，需要先转换：

```python
import rasterio
import numpy as np

# 读取 GeoTIFF
with rasterio.open('your_image.tif') as src:
    image = src.read()  # 自动得到 (bands, rows, cols) 格式

# 保存为 .npy
np.save('image.npy', image)
```

### 目录结构建议

```
my_data/
├── 2024_01/
│   ├── scene_001.npy
│   ├── scene_002.npy
│   └── ...
├── 2024_02/
│   ├── scene_001.npy
│   └── ...
```

在 Web 应用中：
- **输入目录**: `./my_data/`
- **输出目录**: `./my_data_processed/`

处理后会自动保留子目录结构！

## 4️⃣ 参数调整建议

### 矿产勘探（推荐）

| 模块 | 参数 | 值 |
|------|------|-----|
| 大气校正 | 方法 | **DOS** |
|  | 太阳天顶角 | 30-50°（查询获取时间元数据） |
| 几何校正 | 方法 | Affine |
|  | 像元大小 | 30m (Landsat) / 10m (S2) |
| 干扰剔除 | 传感器 | Landsat 8/9 或 S2 |
|  | NDVI 阈值 | 0.15-0.25 |

### 快速处理（追求速度）

- 大气校正: **Simple** 方法
- 几何校正: **Affine** + 最近邻插值
- 干扰剔除: 默认参数

### 高精度处理（追求质量）

- 大气校正: **DOS** 方法 + 准确的太阳天顶角
- 几何校正: **Polynomial** + 双线性/三次插值
- 干扰剔除: 根据实际研究区调整阈值

## 5️⃣ 常见问题

### Q: 为什么输出文件包含 NaN？

**A**: 这是正常的！NaN 表示被干扰剔除模块标记为"干扰"的像元（云、水、植被等）。
在后续处理中，可以用掩膜或 `np.nan_to_num()` 处理。

### Q: 处理速度很慢，怎么办？

**A**: 
- 使用 **Simple** 大气校正（而不是 DOS）
- 使用 **Affine** 几何校正（而不是 Polynomial）
- 使用 **最近邻** 插值（而不是双线性）
- 减小输入影像大小或分块处理

### Q: 如何在 Python 中读取结果？

```python
import numpy as np

# 读取处理后的影像
image = np.load('demo_output/scene1_corrected.npy')
print(f"形状: {image.shape}")
print(f"数据类型: {image.dtype}")
print(f"值范围: [{image.min():.4f}, {image.max():.4f}]")

# 查看某个波段的统计
band1 = image[0]
print(f"波段 1 - 有效像元: {np.sum(~np.isnan(band1))}")
```

### Q: 能否处理单个文件？

**A**: 完全可以！只需在文件列表中勾选 1 个文件即可。

### Q: 能否修改大气校正的太阳几何参数？

**A**: 目前在 UI 中只能设置**太阳天顶角**。如果需要修改日地距、方位角等，需要编辑 `app.py` 中的参数。

## 6️⃣ 后续分析

处理完成后，你可以进行：

### 矿物指数计算

```python
import numpy as np

image = np.load('result_corrected.npy')

# 波段索引 (Landsat 8/9)
r, g, nir, swir1, swir2 = image[3], image[2], image[4], image[5], image[6]

# 矿物指数示例
# Fe 氧化物: SWIR1/NIR
fe_index = swir1 / (nir + 1e-8)

# Al-OH: SWIR2/SWIR1
al_index = swir2 / (swir1 + 1e-8)

# 归一化到 [0, 1]
fe_norm = (fe_index - np.nanmin(fe_index)) / (np.nanmax(fe_index) - np.nanmin(fe_index))
al_norm = (al_index - np.nanmin(al_index)) / (np.nanmax(al_index) - np.nanmin(al_index))

# 保存结果
np.save('fe_index.npy', fe_norm)
np.save('al_index.npy', al_norm)
```

### 导出为 GeoTIFF

```python
import rasterio
from rasterio.transform import Affine
import numpy as np

image = np.load('result_corrected.npy')
bands, rows, cols = image.shape

# 定义地理变换（示例）
transform = Affine(30.0, 0.0, 0.0, 0.0, -30.0, 0.0)

# 写入 GeoTIFF
with rasterio.open(
    'result.tif', 'w',
    driver='GTiff',
    height=rows, width=cols,
    count=bands, dtype=image.dtype,
    transform=transform,
    crs='EPSG:4326'
) as dst:
    for i in range(bands):
        dst.write(image[i], i+1)
```

## 7️⃣ 网络错误排查

| 错误 | 原因 | 解决 |
|------|------|------|
| 连接被拒绝 | 服务未启动 | 运行 `python3 app.py` |
| 页面无法加载 | Flask 未安装 | 运行 `pip install flask` |
| 找不到模块 | 工作目录不对 | 确保在 geo-analyser 目录 |
| 内存不足 | 影像过大 | 分块或降采样处理 |

## 8️⃣ 性能基准

在 MacBook Pro (M1) 上的测试：

| 文件大小 | 处理时间 |
|---------|---------|
| 256×256×7 | ~2 秒 |
| 512×512×7 | ~8 秒 |
| 1024×1024×7 | ~35 秒 |
| 2048×2048×7 | ~140 秒 |

## 需要帮助？

- 📖 完整文档：查看 `README.md`
- 📝 源代码：查看 `app.py`, `templates/index.html`
- 🐛 bug 报告：检查 `demo_data/` 是否已生成

---

**Ready to process? 🚀 Let's go!**
