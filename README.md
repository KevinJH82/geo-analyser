# 🛰 地质蚀变遥感分析平台 (geo-analyser)

面向矿产勘探的 **Flask Web 平台**,聚焦**蚀变异常分析**:从交付影像数据到蚀变成图、成矿预测、形变监测一站式完成。

> 📌 **数据预处理**(大气校正 / 几何校正 / 干扰剔除)已拆分为独立子系统 **`geo-preprocess`**(默认端口 5002),与本系统平级。如需预处理请使用该系统。

## 系统功能

| 模块 | 能力 | 主要接口 |
|------|------|----------|
| 🗂 **交付项目管理** | 加载交付数据、按传感器检测覆盖、上传 ROI(ovkml/geojson)、影像预览 | `/api/projects` `/api/upload_roi` `/api/project_preview` |
| 🔬 **蚀变分析** | 多光谱**波段比值法**与 **Crosta PCA**、高光谱**吸收深度法(band_depth)**;阈值法 `mean_std`/`median_mad`(稳健)/`percentile`;单矿物与批量分析 | `/api/analyze` `/api/analyze_batch` `/api/analyze_preview` |
| 🌿 **丛林模式 + 成矿预测** | 密林区地植物学红边胁迫探测;多证据(蚀变/构造/InSAR形变/冠层)融合成矿有利度;光谱解混与天然出露靶向 | `/api/analyze_batch` (`jungle_mode=true`) |
| 🛰 **多/高光谱传感器** | ASTER、Sentinel-2、Landsat-8(多光谱);EnMAP、PRISMA(高光谱) | — |
| 🧭 **矿床类型推理** | 按矿种推荐蚀变矿物组合,可选 LLM 辅助识别矿床类型 | `/api/deposit_types` `/api/infer_deposit_type` `/api/recommend_targets` |
| 🧱 **构造约束** | 复用 geo-stru 构造解译产物(距断裂邻近度)对靶区加权 | `/api/analyze_batch` (`use_structural_prior=true`) |
| 📈 **InSAR 形变** | 干涉形变分析与时间序列 | `/api/insar/analyze` `/api/insar/timeseries` `/api/insar/stacks` |
| 💾 **结果落盘与下载** | 自动保存 GeoTIFF/PNG/manifest,历史结果浏览与栅格打包下载 | `/api/saved_runs` `/api/download_run_rasters/<run_id>` |

## 启动

```bash
pip install -r requirements.txt
python app.py        # 或 ./run.sh
# 访问 http://127.0.0.1:5001
```

## 关键配置(环境变量)

| 变量 | 默认 | 说明 |
|---|---|---|
| `DELIVERY_ROOT` | (见 delivery_project.py) | 交付项目根目录 |
| `RESULTS_ROOT` | `./results/` | 蚀变结果落盘目录 |
| `DEEPSEEK_API_KEY` | (空) | 矿床类型 LLM 推理(可选) |

## 核心模块

| 文件 | 职责 |
|---|---|
| `alteration_analysis.py` | 波段比值 / Crosta PCA / band_depth / 地植物学胁迫 算法 |
| `alteration_db.py` + `alteration_deposit_db.json` | 矿床类型 ↔ 蚀变矿物 ↔ 波段指标数据库 |
| `prospectivity.py` | 多证据成矿预测融合 |
| `spectral_unmix.py` | 光谱解混 + 天然出露检测 |
| `structural_weighting.py` | geo-stru 构造约束消费 |
| `insar_analysis.py` / `insar_timeseries.py` | InSAR 形变分析 |
| `deposit_type_inference.py` | LLM 辅助矿床类型识别 |
| `alteration_store.py` | 结果持久化(GeoTIFF/PNG/manifest) |
| `delivery_project.py` | 交付项目 / 多传感器数据装载 |

> 纯光谱指数函数(NDVI/NDRE/BSI 等)统一来自 `commons/spectral_indices.py`,与 geo-preprocess 同源共享。

## 相关系统(同级)
- **geo-preprocess** — 数据预处理(大气/几何/干扰剔除)
- **geo-downloader** — 遥感数据下载
- **geo-stru** — 构造解译(经 `commons/structural_broker` 消费)
- **geo-insar** — InSAR 形变(经 `commons/insar_broker` 消费)
