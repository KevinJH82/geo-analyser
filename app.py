"""
遥感数据预处理 Web 应用
集成：大气校正 → 几何校正 → 干扰剔除
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import io
import base64
from pathlib import Path
from flask import Flask, render_template, request, jsonify, Response, send_file
from datetime import datetime
import traceback

# 从同级 .env 加载环境变量(如 DEEPSEEK_API_KEY),让密钥持久化、无需每次手动 export。
# 必须在导入会读取 env 的项目模块之前调用。优先 python-dotenv,缺失则用无依赖手动解析
# (系统 python 常无 dotenv,此前会静默忽略 .env)。
_envf = Path(__file__).parent / ".env"
try:
    from dotenv import load_dotenv
    load_dotenv(_envf)
except ImportError:
    if _envf.exists():
        for _line in _envf.read_text(encoding="utf-8").splitlines():
            _line = _line.strip()
            if _line and not _line.startswith("#") and "=" in _line:
                _k, _v = _line.split("=", 1)
                os.environ.setdefault(_k.strip(), _v.strip())

from alteration_analysis import (
    analyze_alteration, get_supported_minerals, MINERAL_CATALOG,
    analyze_single, analyze_batch, build_roi_mask, AlterationResult, BatchResult,
)
from alteration_db import (
    list_commodities, list_deposit_types,
    get_recommended_targets, get_targets_multi_sensor, get_deposit_type_meta,
)
from deposit_type_inference import infer_deposit_types_detailed
from structural_deposit import structural_deposit_candidates
from delivery_project import (
    DELIVERY_ROOT, list_projects, open_project_from_upload,
    resolve_project_dir, load_sensor_data, list_available_sensors,
    bbox_from_geojson, check_sensor_coverage,
    load_enmap_data, ENMAP_KEY, get_enmap_metadata_path,
    load_prisma_data, PRISMA_KEY, get_prisma_metadata_paths,
)
from basemap_provider import fetch_satellite_basemap
from alteration_store import (
    save_batch_run, list_runs, load_manifest, resolve_file, results_root,
    make_run_rasters_zip,
)

# matplotlib 使用无 GUI 后端
matplotlib.use('Agg')

app = Flask(__name__, template_folder='templates', static_folder='static')
# ── 内部鉴权:拒绝绕过 BFF 的直连(PORTAL_INTERNAL_KEY 配置后生效) ──
try:
    import sys as _ia_sys
    if '/opt/deepexplor-services' not in _ia_sys.path:
        _ia_sys.path.insert(0, '/opt/deepexplor-services')
    from commons.internal_auth import init_internal_auth as _init_internal_auth
    _init_internal_auth(app)
except Exception as _ia_e:
    print(f'[internal_auth] 跳过接入: {_ia_e}')
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB

# ─────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────

SUPPORTED_EXTENSIONS = {".tif", ".tiff", ".npy"}

# 单波段文件名模式：B1.tif, B2.tif, B3N.tif, B3B.tif, band1.tif 等
# group(1)=数字部分, group(2)=可选字母后缀（如 N/B）
import re
from typing import Optional, Dict, Any
BAND_FILE_PATTERN = re.compile(r'^[Bb](?:and)?(\d+)([A-Za-z]*)\.(tif|tiff)$', re.IGNORECASE)


def _band_sort_key(name: str):
    """返回 (数字, 字母后缀) 用于波段排序，如 B3N → (3, 'N')"""
    m = BAND_FILE_PATTERN.match(name)
    if not m:
        return (9999, name)
    return (int(m.group(1)), m.group(2).upper())


def get_band_files(directory: Path) -> list:
    """
    检测目录中的波段文件（B1.tif, B2.tif, B3N.tif ...），
    返回按波段号排序的文件路径列表。
    """
    band_files = [f for f in directory.iterdir() if BAND_FILE_PATTERN.match(f.name)]
    if not band_files:
        return []
    return sorted(band_files, key=lambda f: _band_sort_key(f.name))


def read_image(image_path: str) -> tuple:
    """
    读取影像，返回 (image_array, profile)。
    - 若 image_path 是目录：按 B1.tif, B2.tif ... 顺序合并为多波段
    - 若是单 .tif/.tiff 文件：直接读取
    - 若是 .npy：np.load
    """
    import rasterio
    p = Path(image_path)

    if p.is_dir():
        band_files = get_band_files(p)
        if not band_files:
            raise ValueError(f"目录 {p} 中未找到波段文件（B1.tif, B2.tif ...）")

        # 按分辨率分组，取数量最多的分辨率组，组内裁剪到最小公共尺寸
        res_groups = {}
        for bf in band_files:
            with rasterio.open(bf) as src:
                res = round(abs(src.res[0]))
            res_groups.setdefault(res, []).append(bf)
        # 选数量最多的分辨率组（通常是 10m 或 30m 光学波段）
        dominant_res = max(res_groups, key=lambda r: len(res_groups[r]))
        band_files = res_groups[dominant_res]

        # 读出各波段尺寸，裁剪到最小公共 rows/cols
        band_shapes = {}
        for bf in band_files:
            with rasterio.open(bf) as src:
                band_shapes[bf] = (src.height, src.width)
        min_rows = min(s[0] for s in band_shapes.values())
        min_cols = min(s[1] for s in band_shapes.values())

        bands = []
        profile = None
        for bf in band_files:
            with rasterio.open(bf) as src:
                bands.append(src.read(1)[:min_rows, :min_cols].astype(np.float32))
                if profile is None:
                    profile = src.profile.copy()
        image = np.stack(bands, axis=0)  # (bands, rows, cols)
        profile.update(count=len(bands), height=min_rows, width=min_cols)
        return image, profile

    suffix = p.suffix.lower()
    if suffix in (".tif", ".tiff"):
        with rasterio.open(image_path) as src:
            image = src.read().astype(np.float32)
            profile = src.profile.copy()
        return image, profile

    return np.load(image_path, allow_pickle=True), None


# ─────────────────────────────────────────────
# Flask 路由
# ─────────────────────────────────────────────

@app.route('/')
def index():
    """主页"""
    from flask import make_response
    resp = make_response(render_template('index.html'))
    resp.headers['Cache-Control'] = 'no-store'
    return resp


@app.route('/docs/alteration-reference')
def docs_alteration_reference():
    """在线浏览蚀变遥感对照手册 — md 文件用 marked.js 渲染。"""
    import json as _json
    md_path = Path(__file__).parent / 'docs' / 'alteration_remote_sensing_reference.md'
    try:
        md_text = md_path.read_text(encoding='utf-8')
    except Exception as e:
        return Response(f"文档读取失败: {e}", status=500, mimetype='text/plain; charset=utf-8')
    md_js = _json.dumps(md_text)
    html = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<title>蚀变遥感对照手册 · DeepExplor</title>
<script src="https://cdn.jsdelivr.net/npm/marked/marked.min.js"></script>
<style>
  * { box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, "PingFang SC", "Microsoft YaHei", sans-serif; color: #1e293b; line-height: 1.7; max-width: 1100px; margin: 0 auto; padding: 32px 40px; background: #fafbfc; }
  h1 { font-size: 28px; margin: 0 0 16px; padding-bottom: 10px; border-bottom: 3px solid #059669; }
  h2 { font-size: 20px; margin: 32px 0 12px; padding-bottom: 6px; border-bottom: 2px solid #059669; color: #059669; }
  h3 { font-size: 16px; margin: 20px 0 8px; color: #334155; }
  h4 { font-size: 14px; margin: 14px 0 6px; color: #475569; }
  table { width: 100%; border-collapse: collapse; font-size: 12px; margin: 10px 0; background: white; }
  th, td { padding: 7px 10px; border: 1px solid #e2e8f0; text-align: left; vertical-align: top; }
  th { background: #f1f5f9; font-weight: 600; }
  code { background: #f1f5f9; padding: 1px 6px; border-radius: 3px; font-family: "Menlo", "Consolas", monospace; font-size: 0.9em; color: #be123c; }
  pre code { display: block; padding: 12px 16px; background: #1e293b; color: #f1f5f9; border-radius: 6px; overflow-x: auto; }
  blockquote { border-left: 3px solid #94a3b8; padding-left: 14px; color: #475569; margin: 12px 0; background: #f8fafc; padding-top: 8px; padding-bottom: 8px; border-radius: 0 4px 4px 0; }
  hr { border: none; border-top: 1px solid #cbd5e1; margin: 24px 0; }
  a { color: #059669; text-decoration: none; }
  a:hover { text-decoration: underline; }
  ul, ol { padding-left: 26px; }
  li { margin: 3px 0; }
  .toolbar { position: sticky; top: 0; background: rgba(250,251,252,0.95); backdrop-filter: blur(6px); padding: 12px 0; margin: -12px 0 8px; border-bottom: 1px solid #e2e8f0; z-index: 10; }
  .toolbar a { font-size: 12px; padding: 6px 12px; background: white; border: 1px solid #cbd5e1; border-radius: 6px; margin-right: 8px; }
  .toolbar a:hover { background: #f1f5f9; text-decoration: none; }
</style></head>
<body>
<div class="toolbar">
  <a href="javascript:window.print();">🖨 打印 / 存 PDF</a>
  <a href="/docs/alteration-reference.md" download="alteration_remote_sensing_reference.md">📥 下载 Markdown 原文</a>
  <a href="/" target="_blank">↗ 回主页</a>
</div>
<div id="content"></div>
<script>
  const md = __MD_TEXT__;
  document.getElementById('content').innerHTML = marked.parse(md, { gfm: true, breaks: false });
</script>
</body></html>"""
    html = html.replace("__MD_TEXT__", md_js)
    return Response(html, mimetype='text/html; charset=utf-8')


@app.route('/docs/alteration-reference.md')
def docs_alteration_reference_md():
    """下载 markdown 原文。"""
    md_path = Path(__file__).parent / 'docs' / 'alteration_remote_sensing_reference.md'
    try:
        md_text = md_path.read_text(encoding='utf-8')
    except Exception as e:
        return Response(f"文档读取失败: {e}", status=500, mimetype='text/plain; charset=utf-8')
    return Response(md_text, mimetype='text/markdown; charset=utf-8')


# ─────────────────────────────────────────────
# 蚀变分析路由
# ─────────────────────────────────────────────

@app.route('/api/mineral_list', methods=['POST'])
def api_mineral_list():
    """
    返回指定传感器支持的矿物列表。
    请求体: {"sensor": "Landsat8/9_L2"} 或 {"sensor": "ASTER"}
    """
    data = request.json or {}
    sensor_raw = data.get("sensor", "Landsat8/9")
    sensor = sensor_raw.replace("_L2", "")   # 去掉 L2 后缀统一匹配
    minerals = get_supported_minerals(sensor)
    return jsonify({"minerals": minerals, "sensor": sensor})


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    """
    对指定影像执行蚀变分析，返回结果统计 + 叠加预览图（base64 PNG）。

    请求体:
    {
        "file_path": "/path/to/image.tif",   # 文件或波段目录
        "sensor": "Landsat8/9_L2",
        "mineral": "iron_oxide",
        "threshold_method": "mean_std",       # "mean_std" | "median_mad" | "percentile"
        "k": 1.5,
        "colormap": "hot"                     # matplotlib colormap 名
    }
    """
    data = request.json or {}
    file_path = data.get("file_path", "")
    sensor_raw = data.get("sensor", "Landsat8/9")
    sensor = sensor_raw.replace("_L2", "")
    mineral = data.get("mineral", "iron_oxide")
    threshold_method = data.get("threshold_method", "mean_std")
    k = float(data.get("k", 1.5))
    colormap = data.get("colormap", "hot")
    bg_image_b64 = data.get("bg_image", None)  # 用户上传的底图（base64）

    if not file_path or (not os.path.isfile(file_path) and not os.path.isdir(file_path)):
        return jsonify({"error": "文件不存在"}), 400

    try:
        image, profile = read_image(file_path)

        # L2 整数格式归一化
        if sensor_raw.endswith("_L2") and image.max() > 10.0:
            image = image / 10000.0

        result = analyze_alteration(image, sensor, mineral, threshold_method, k)

        # 生成叠加可视化图（base64 PNG）
        if bg_image_b64:
            overlay_b64 = _render_overlay_on_bg(bg_image_b64, result, colormap)
        else:
            overlay_b64 = _render_alteration_overlay(image, result, colormap)

        return jsonify({
            "mineral": result.mineral,
            "label": result.label,
            "sensor": result.sensor,
            "anomaly_ratio": round(result.anomaly_ratio * 100, 2),  # 百分比
            "threshold": float(result.threshold) if np.isfinite(result.threshold) else None,
            "warning": result.warning,
            "overlay": f"data:image/png;base64,{overlay_b64}" if overlay_b64 else None,
        })

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": f"分析失败: {str(e)}"}), 500


# ─────────────────────────────────────────────
# 新蚀变 API (v2): 矿床类型驱动 + 双方法批量
# ─────────────────────────────────────────────

# ─────────────────────────────────────────────
# 蚀变 v3: 交付数据项目驱动
# ─────────────────────────────────────────────

UPLOAD_TMP = Path("/tmp/geo_analyser_roi_uploads")
UPLOAD_TMP.mkdir(parents=True, exist_ok=True)


@app.route('/api/projects', methods=['GET'])
def api_projects():
    """列出交付目录下所有项目(供前端兜底浏览)。"""
    return jsonify({
        "delivery_root": str(DELIVERY_ROOT),
        "exists":        DELIVERY_ROOT.is_dir(),
        "projects":      list_projects(),
    })


@app.route('/api/upload_roi', methods=['POST'])
def api_upload_roi():
    """
    用户上传一个 ROI 文件(ovkml/kml/geojson)。
    按文件主名称在交付目录中定位项目,解析 ROI,扫描可用传感器。
    返回 {project_name, project_dir, roi_geojson, bbox, available_sensors}
    或 {error, ...}
    """
    if "file" not in request.files:
        return jsonify({"error": "缺少上传文件(form field 'file')"}), 400
    f = request.files["file"]
    if not f.filename:
        return jsonify({"error": "文件名为空"}), 400

    saved = UPLOAD_TMP / f.filename
    f.save(saved)

    info = open_project_from_upload(saved)
    if info is None or "error" in info:
        msg = info.get("error", "未知错误") if info else "解析失败"
        return jsonify({"error": msg,
                        "delivery_root": str(DELIVERY_ROOT)}), 404
    return jsonify(info)


@app.route('/api/project_preview', methods=['POST'])
def api_project_preview():
    """
    返回 ROI 预览图 PNG。优先用在线卫星瓦片底图 + ROI 多边形描边;
    在线失败时 fallback 到该传感器的真彩色 RGB。
    body: {project_name, sensor_key, roi_geojson}
    """
    data = request.json or {}
    project_name = data.get("project_name", "")
    sensor_key   = data.get("sensor_key", "")
    roi_geojson  = data.get("roi_geojson")
    if not project_name or not sensor_key:
        return jsonify({"error": "缺少 project_name 或 sensor_key"}), 400

    project_dir = DELIVERY_ROOT / project_name
    if not project_dir.is_dir():
        return jsonify({"error": f"项目目录不存在: {project_name}"}), 404

    from PIL import Image as PILImage
    # 智能选底图: 项目内 Sentinel-2 → ASTER → Esri
    basemap_used = None
    pil = None
    if roi_geojson:
        # 计算项目可用且覆盖 ROI 的传感器
        sensors_info = list_available_sensors(project_dir)
        avail = {s["key"] for s in sensors_info}
        from delivery_project import check_sensor_coverage as _ck
        roi_bbox = bbox_from_geojson(roi_geojson)
        usable = set()
        if roi_bbox:
            for sk in avail:
                if _ck(project_dir, sk, roi_bbox).get("covers"):
                    usable.add(sk)
        bm = _make_basemap_view(roi_geojson, target_max_px=2400,
                                project_dir=project_dir,
                                available_sensors=usable)
        if bm is not None:
            # 在底图上叠加 ROI 描边(用 PIL)
            img = bm["image"].copy()
            try:
                from PIL import ImageDraw
                pil = PILImage.fromarray(img, mode="RGB").convert("RGBA")
                draw = ImageDraw.Draw(pil)
                from rasterio.transform import rowcol
                geom = roi_geojson["geometry"] if "geometry" in roi_geojson else roi_geojson
                if geom.get("type") == "Polygon":
                    for ring in geom["coordinates"]:
                        pts = []
                        for lon, lat, *_ in ring:
                            row, col = rowcol(bm["transform"], lon, lat)
                            pts.append((col, row))
                        if len(pts) >= 2:
                            draw.line(pts, fill=(16, 185, 129, 255), width=3)
                pil = pil.convert("RGB")
                basemap_used = bm.get("source", "satellite")
            except Exception as e:
                app.logger.warning(f"ROI 描边失败: {e}")
                pil = PILImage.fromarray(img, mode="RGB")
                basemap_used = bm.get("source", "satellite")

    if pil is None:
        # fallback: 传感器自家真彩色
        try:
            if sensor_key in (ENMAP_KEY, PRISMA_KEY):
                # 高光谱(EnMAP/PRISMA):按波长取近似真彩色波段 (R~0.66 G~0.55 B~0.47 µm)
                _hs_loader = load_enmap_data if sensor_key == ENMAP_KEY else load_prisma_data
                cube, wl, _ = _hs_loader(project_dir)
                def _nearest(target):
                    return cube[int(np.argmin(np.abs(wl - target)))]
                rgb = np.dstack([_nearest(0.66), _nearest(0.55), _nearest(0.47)])
                # 2–98% 拉伸增强对比
                finite = rgb[np.isfinite(rgb)]
                if finite.size:
                    lo, hi = np.nanpercentile(rgb, 2), np.nanpercentile(rgb, 98)
                    if hi > lo:
                        rgb = (rgb - lo) / (hi - lo)
                rgb = np.nan_to_num(np.clip(rgb, 0, 1))
            else:
                image, bn_map, profile = load_sensor_data(project_dir, sensor_key)
                if image.max() > 10.0:
                    image = image / 10000.0
                rgb = _make_rgb_base_dyn(image, sensor_key, bn_map)
            rgb_uint8 = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
            pil = PILImage.fromarray(rgb_uint8, mode="RGB")
            basemap_used = f"sensor_rgb:{sensor_key}"
        except Exception as e:
            return jsonify({"error": str(e)}), 500

    max_side = 1024
    W, H = pil.size
    scale = min(max_side / max(H, W), 1.0)
    if scale < 1.0:
        pil = pil.resize((int(W * scale), int(H * scale)), PILImage.BILINEAR)

    buf = io.BytesIO()
    pil.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("utf-8")

    return jsonify({
        "preview":      f"data:image/png;base64,{b64}",
        "width":        W,
        "height":       H,
        "preview_w":    pil.size[0],
        "preview_h":    pil.size[1],
        "sensor_key":   sensor_key,
        "basemap":      basemap_used,
    })


@app.route('/api/raw_preview', methods=['POST'])
def api_raw_preview():
    """
    返回影像的 1:1 RGB 预览 PNG (用 PIL 渲染,无 matplotlib margin),
    便于前端在画布上以精确像素坐标绘制 ROI。
    body: {file_path, sensor}
    返回: {preview, width, height, has_georef}
    """
    data = request.json or {}
    file_path = data.get("file_path", "")
    sensor_raw = data.get("sensor", "Landsat8/9")
    sensor = sensor_raw.replace("_L2", "")
    if not file_path or (not os.path.isfile(file_path) and not os.path.isdir(file_path)):
        return jsonify({"error": "文件不存在"}), 400
    try:
        image, profile = read_image(file_path)
        if sensor_raw.endswith("_L2") and image.max() > 10.0:
            image = image / 10000.0

        rgb = _make_rgb_base(image, sensor)
        rgb_uint8 = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)

        from PIL import Image as PILImage
        pil = PILImage.fromarray(rgb_uint8, mode="RGB")
        # 限制最大边到 1024 像素,加速传输
        max_side = 1024
        H, W = rgb_uint8.shape[:2]
        scale = min(max_side / max(H, W), 1.0)
        if scale < 1.0:
            pil = pil.resize((int(W * scale), int(H * scale)), PILImage.BILINEAR)

        buf = io.BytesIO()
        pil.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")

        has_georef = bool(profile and profile.get("transform") and profile.get("crs"))
        return jsonify({
            "preview":     f"data:image/png;base64,{b64}",
            "width":       W,
            "height":      H,
            "preview_w":   pil.size[0],
            "preview_h":   pil.size[1],
            "has_georef":  has_georef,
        })
    except Exception as e:
        app.logger.error(f"raw_preview 失败: {e}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@app.route('/api/commodity_list', methods=['GET'])
def api_commodity_list():
    """27 个矿种总览。"""
    return jsonify({"commodities": list_commodities()})


@app.route('/api/deposit_types', methods=['POST'])
def api_deposit_types():
    """指定矿种下的矿床类型列表。"""
    data = request.json or {}
    commodity = data.get("commodity", "")
    if not commodity:
        return jsonify({"error": "缺少 commodity"}), 400
    return jsonify({
        "commodity":      commodity,
        "deposit_types":  list_deposit_types(commodity),
    })


@app.route('/api/recommend_targets', methods=['POST'])
def api_recommend_targets():
    """
    根据矿床类型 + (可选)项目可用传感器,返回推荐蚀变矿物列表的多传感器视图。
    body: {deposit_type, available_sensors?: ["ASTER","Sentinel2"]}

    每个 target 标注:
      - preferred_sensor: 库推荐的首选传感器
      - effective_sensor: 综合项目可用情况后实际会跑的传感器(可能 fallback,可能为 None)
      - per_sensor: 每个传感器的 ratio/pca 可用性
      - data_present: 该 mineral 的 effective_sensor 是否在 available_sensors 中
    """
    data = request.json or {}
    deposit_type = data.get("deposit_type", "")
    available = set(data.get("available_sensors") or [])
    if not deposit_type:
        return jsonify({"error": "缺少 deposit_type"}), 400

    multi = get_targets_multi_sensor(deposit_type)
    meta = get_deposit_type_meta(deposit_type) or {}

    targets = []
    for t in multi:
        preferred = t["preferred_sensor"]
        # effective_sensor: preferred 在 available 内则用 preferred;
        # 否则按 ASTER>Sentinel2>Landsat8 顺序找第一个 available 且至少一种方法可用的
        effective = None
        if preferred and (not available or preferred in available):
            effective = preferred
        if effective is None and available:
            for s in ("ASTER", "Sentinel2", "Landsat8"):
                if s in available:
                    ps = t["per_sensor"].get(s, {})
                    if ps.get("ratio_available") or ps.get("pca_available"):
                        effective = s
                        break
        # 高光谱兜底:多光谱无可用项时,若项目有 EnMAP/PRISMA 且该目标支持吸收深度法(band_depth),
        # 用高光谱。这样"高光谱烃指数"等多光谱无对应波段的目标在有高光谱数据时也能被推荐、可勾选。
        if effective is None and available:
            for s in (ENMAP_KEY, PRISMA_KEY):
                if s in available and t["per_sensor"].get(s, {}).get("band_depth_available"):
                    effective = s
                    break

        targets.append({
            "mineral":          t["mineral"],
            "zone":             t["zone"],
            "priority":         t["priority"],
            "anomaly_type":     t["anomaly_type"],
            "preferred_sensor": preferred,
            "effective_sensor": effective,
            "data_present":     bool(effective and effective in available) if available else None,
            "per_sensor": {
                s: {
                    "ratio_expr":           ps["ratio_expr"],
                    "ratio_available":      ps["ratio_available"],
                    "pca_available":        ps["pca_available"],
                    "pca_bands":            (ps["pca_spec"] or {}).get("input_bands", []),
                    "band_depth_available": ps.get("band_depth_available", False),
                }
                for s, ps in t["per_sensor"].items()
            },
        })

    return jsonify({
        "deposit_type":      deposit_type,
        "available_sensors": sorted(available) if available else [],
        "meta":              meta,
        "targets":           targets,
    })


def _pixel_polygon_to_geojson(file_path: str, pixel_polygon: list) -> Optional[Dict[str, Any]]:
    """
    把像素坐标多边形 [[x,y],...] 用影像的 GeoTransform 转成 EPSG:4326 GeoJSON。
    无 transform 或转换失败返回 None。
    """
    try:
        _, profile = read_image(file_path)
        if not profile:
            return None
        transform = profile.get("transform")
        crs = profile.get("crs")
        if transform is None:
            return None

        # 像素 → 源坐标系
        src_coords = [transform * (float(x), float(y)) for x, y in pixel_polygon]
        # 闭合环
        if src_coords[0] != src_coords[-1]:
            src_coords.append(src_coords[0])

        # 若 CRS 已是 EPSG:4326,直接用;否则 reproject
        if crs is None:
            return None
        if crs.to_epsg() == 4326:
            ring = [[lon, lat] for lon, lat in src_coords]
        else:
            from rasterio.warp import transform as warp_transform
            xs = [c[0] for c in src_coords]
            ys = [c[1] for c in src_coords]
            lons, lats = warp_transform(crs, "EPSG:4326", xs, ys)
            ring = [[lon, lat] for lon, lat in zip(lons, lats)]

        return {"type": "Polygon", "coordinates": [ring]}
    except Exception as e:
        app.logger.warning(f"_pixel_polygon_to_geojson 失败: {e}")
        return None


@app.route('/api/structural_deposit_type', methods=['POST'])
def api_structural_deposit_type():
    """
    最高优先级来源:从 geo-stru 构造解译产出获取矿床类型候选(翻译成蚀变库 type_name)。

    body: {roi_geojson: {...}, top_k: 3}
    返回: {candidates: [{deposit_type, confidence, evidence, source}],
           degraded: bool, hint: str, source: "geo-stru", aoi_name, run_id}
    无 geo-stru 产物 / 映射不到时 degraded=true,前端据此回退 LLM 自动识别。
    任何异常静默降级,绝不阻断前端级联。
    """
    data = request.json or {}
    top_k = int(data.get("top_k", 3))
    roi = data.get("roi_geojson")
    if not roi:
        return jsonify({
            "candidates": [], "degraded": True, "source": "geo-stru",
            "hint": "缺少 roi_geojson",
        })
    try:
        result = structural_deposit_candidates(roi, top_k=top_k)
        candidates = result["candidates"]
        status = result.get("status", "no_product")
        return jsonify({
            "candidates": candidates,
            "status":     status,                       # ok | no_product | non_alteration
            "degraded":   len(candidates) == 0,
            # non_alteration:geo-stru 已权威判定为油气/非蚀变靶区,前端不应回退 LLM 改判
            "applicable": status != "non_alteration",
            "source":     "geo-stru",
            "aoi_name":   result.get("aoi_name"),
            "run_id":     result.get("run_id"),
            "primary_model": result.get("primary_model"),
            "primary_confidence": result.get("primary_confidence"),
            "hint":       result.get("reason") or "",
        })
    except Exception as e:
        app.logger.error(f"structural_deposit_type 失败: {e}", exc_info=True)
        return jsonify({
            "candidates": [], "status": "no_product", "degraded": True,
            "applicable": True, "source": "geo-stru",
            "hint": f"geo-stru 构造推理失败: {e}",
        })


@app.route('/api/infer_deposit_type', methods=['POST'])
def api_infer_deposit_type():
    """
    根据 ROI 自动识别 top-3 候选矿床类型。

    body 接受两种形式:
      {roi_geojson: {...}, top_k: 3}          - 直接传 GeoJSON
      {file_path, pixel_polygon, top_k: 3}    - 传像素 ROI + 文件路径,后端用 GeoTransform 转 GeoJSON

    返回: {candidates: [{deposit_type, confidence, evidence, source}],
           degraded: bool, hint: str, georeferenced: bool}
    无地理参考时 degraded=true,前端进入手动模式。
    """
    data = request.json or {}
    top_k = int(data.get("top_k", 3))

    roi = data.get("roi_geojson")
    if not roi:
        pixel_polygon = data.get("pixel_polygon")
        file_path = data.get("file_path", "")
        if not pixel_polygon or not file_path:
            return jsonify({"error": "缺少 roi_geojson 或 (file_path + pixel_polygon)"}), 400
        roi = _pixel_polygon_to_geojson(file_path, pixel_polygon)
        if roi is None:
            return jsonify({
                "candidates": [], "degraded": True, "georeferenced": False,
                "hint": "影像无地理参考,无法自动识别。请手动选择矿种和矿床类型。",
            })

    try:
        result = infer_deposit_types_detailed(roi, top_k=top_k)
        candidates = result["candidates"]
        return jsonify({
            "candidates": candidates,
            "degraded":      len(candidates) == 0,
            "georeferenced": True,
            "hint":          (result.get("reason") or "自动识别未返回结果,请手动选择") if not candidates else "",
        })
    except Exception as e:
        app.logger.error(f"infer_deposit_type 失败: {e}", exc_info=True)
        return jsonify({
            "candidates": [], "degraded": True, "georeferenced": True,
            "hint": str(e),
        })


def _apply_structural_prior(score, shape, transform, crs, roi_geojson, weight=0.15):
    """
    用 geo-stru 构造产物(距断裂邻近度)对蚀变综合得分加权:成矿受构造控制,
    近断裂处上调、远断裂轻度下调,把分散异常约束为沿构造展布的靶区。

    复用 geo-analyser/structural_weighting + commons/structural_broker(bbox 相交发现)。
    任何缺失/失败都原样返回 score——可选增强,绝不破坏既有分析。
    """
    try:
        import sys
        if '/opt/deepexplor-services' not in sys.path:
            sys.path.insert(0, '/opt/deepexplor-services')
        import structural_weighting as _sw
        if not roi_geojson or transform is None or crs is None:
            return score
        bbox = bbox_from_geojson(roi_geojson)
        if not bbox:
            return score
        layers = _sw.load_structural_layers(bbox, shape, transform, crs)
        if not layers or layers.get('distance') is None:
            return score
        prox = _sw.proximity_from_distance(layers['distance'])
        out = _sw.apply_structural_weighting(score, prox, weight=weight)
        app.logger.info(f"已叠加 geo-stru 构造控矿先验 (AOI={layers.get('aoi_name')}, weight={weight})")
        return out
    except Exception as e:
        app.logger.warning(f"构造控矿先验叠加失败,忽略: {e}")
        return score


def _or_masks(method_masks: Dict[str, np.ndarray]):
    """把一个矿物在多方法下的异常掩膜取并集(任一方法判异常即为异常)。"""
    arrs = [v for v in (method_masks or {}).values() if v is not None]
    if not arrs:
        return None
    u = arrs[0].astype(bool).copy()
    for x in arrs[1:]:
        u |= x.astype(bool)
    return u


def _structural_assoc_for_sensor(bbox, shape, transform, crs, masks_by_mineral, buffer_m=300.0):
    """对某传感器格网载入'距断裂距离'层,逐矿物算异常-构造关联度。返回 {mineral: assoc}。"""
    out: Dict[str, Any] = {}
    if not bbox or transform is None or crs is None or not masks_by_mineral:
        return out
    try:
        import structural_weighting as _sw
        layers = _sw.load_structural_layers(bbox, shape, transform, crs)
        dist = (layers or {}).get("distance")
        if dist is None:
            return out
        for mineral, mask in masks_by_mineral.items():
            a = _sw.lineament_association(mask, dist, buffer_m=buffer_m)
            if a is not None:
                out[mineral] = a
    except Exception as e:
        app.logger.warning(f"构造关联度计算失败,忽略: {e}")
    return out


def _compute_prospectivity(score_map, veg_stress_index, shape, transform, crs,
                           roi_bbox, roi_mask, weights=None):
    """
    丛林模式 capstone:收集多证据层融合为成矿有利度。

      E1 alteration   = 蚀变综合得分 score_map
      E2 veg_stress   = 地植物学胁迫指数(红边)
      E3 structure    = 构造邻近度(geo-stru 距断裂距离 → 指数衰减)
      E4 deformation  = InSAR 形变证据(geo-insar deformation_evidence,重投影到本格网)

    任何证据层缺失(数据不存在/重投影失败)静默忽略,prospectivity.fuse_evidence 自适应再归一权重。
    返回 prospectivity.ProspectivityResult,失败返回 None。
    """
    try:
        import prospectivity as _pr
    except Exception as e:
        app.logger.warning(f"prospectivity 模块不可用: {e}")
        return None

    layers = {"alteration": score_map, "veg_stress": veg_stress_index}

    # E3 构造邻近度
    try:
        import structural_weighting as _sw
        if roi_bbox and transform is not None and crs is not None:
            sl = _sw.load_structural_layers(roi_bbox, shape, transform, crs)
            if sl and sl.get("distance") is not None:
                layers["structure"] = _sw.proximity_from_distance(sl["distance"])
    except Exception as e:
        app.logger.info(f"prospectivity 构造证据层缺失,忽略: {e}")

    # E4 InSAR 形变证据
    try:
        from commons.insar_broker import find_insar_for_bbox, get_product_path
        import structural_weighting as _sw
        if roi_bbox and transform is not None and crs is not None:
            matches = find_insar_for_bbox(roi_bbox)
            if matches:
                dpath = get_product_path(matches[0], "deformation_evidence")
                if dpath:
                    layers["deformation"] = _sw._reproject_to(dpath, shape, transform, crs)
    except Exception as e:
        app.logger.info(f"prospectivity InSAR 形变证据层缺失,忽略: {e}")

    return _pr.fuse_evidence(layers, weights=weights, roi_mask=roi_mask)


def _png_to_datauri(path) -> Optional[str]:
    """把本地 PNG 文件读成 data URI(供前端报告内嵌走向玫瑰图)。失败返回 None。"""
    try:
        if not path:
            return None
        with open(path, "rb") as f:
            return "data:image/png;base64," + base64.b64encode(f.read()).decode("utf-8")
    except Exception:
        return None


def _run_hyperspectral_analysis(
    sensor_key, cube, wavelengths_um, profile, aux_paths,
    targets_used, roi_geojson, basemap_view,
    threshold_method, k, colormap,
    per_mineral_results, composites_per_sensor, save_sensors,
    use_structural_prior=False,
    struct_lineaments=None, struct_intersections=None,
    roi_bbox=None, struct_assoc_out=None,
):
    """
    对一个高光谱传感器(EnMAP / PRISMA)跑吸收深度(band_depth)分析,
    把逐矿物结果与综合判断写入 per_mineral_results / composites_per_sensor / save_sensors。

    与多光谱的 ratio/pca 并列,作为额外的传感器列。返回 (high_confidence_pixels, roi_size)。
    """
    targets = [t for t in targets_used
               if t["per_sensor"].get(sensor_key, {}).get("band_depth_available")]
    if not targets:
        return 0, 0

    Hh, Ww = cube.shape[1], cube.shape[2]
    transform = profile.get("transform")
    crs = profile.get("crs")
    roi_mask = build_roi_mask((Hh, Ww), roi_geojson=roi_geojson,
                              transform=transform, image_crs=crs)
    roi_size = int(roi_mask.sum()) if roi_mask is not None else Hh * Ww

    # 源数据的辅助文件(EnMAP METADATA.XML / PRISMA *.hdr)随结果落盘,供按波长解译
    save_sensors[sensor_key] = {
        "profile": profile, "shape": [Hh, Ww], "roi_pixels": roi_size,
        "results": [], "composite": None,
        "aux_files": [str(p) for p in (aux_paths or []) if p],
    }
    masks: Dict[str, np.ndarray] = {}

    for t in targets:
        mineral = t["mineral"]
        entry = per_mineral_results.setdefault(mineral, {
            "mineral":          mineral,
            "zone":             t["zone"],
            "priority":         t["priority"],
            "anomaly_type":     t["anomaly_type"],
            "preferred_sensor": t["preferred_sensor"],
            "effective_sensor": t.get("effective_sensor"),
            "data_status":      "ok" if t.get("effective_sensor") else "hyperspectral_only",
            "notes":            list(t.get("notes") or []),
            "results":          {},
        })
        sensor_results = entry["results"].setdefault(sensor_key, {})
        feat = t["per_sensor"][sensor_key].get("enmap_feature") or {}
        target_spec = {
            "mineral":              mineral,
            "enmap_feature":        feat,
            "band_depth_available": True,
        }
        try:
            res = analyze_single(
                cube, sensor_key, target_spec, "band_depth",
                roi_mask=roi_mask, threshold_method=threshold_method,
                k=k, wavelengths_um=wavelengths_um,
            )
            if basemap_view is not None:
                cell_b64 = _render_cell_on_basemap(
                    basemap_view, res.anomaly_mask,
                    transform, crs, roi_geojson, res, sensor_key, colormap,
                )
            else:
                cell_b64 = None
            sensor_results["band_depth"] = {
                "anomaly_ratio": round(res.anomaly_ratio * 100, 2),
                "threshold":     None if not np.isfinite(res.threshold) else float(res.threshold),
                "warning":       res.warning,
                "overlay":       f"data:image/png;base64,{cell_b64}" if cell_b64 else None,
                "ratio_expr":    res.ratio_expr,
                "sign":          res.sign,
                "feature_um":    feat.get("feature_um"),
            }
            masks[mineral] = res.anomaly_mask
            save_sensors[sensor_key]["results"].append({
                "mineral":          mineral,
                "zone":             t["zone"],
                "priority":         t["priority"],
                "anomaly_type":     t["anomaly_type"],
                "effective_sensor": sensor_key,
                "data_status":      entry["data_status"],
                "method":           "band_depth",
                "anomaly_ratio":    round(res.anomaly_ratio * 100, 2),
                "threshold":        None if not np.isfinite(res.threshold) else float(res.threshold),
                "ratio_expr":       res.ratio_expr,
                "pc_used":          None,
                "sign":             res.sign,
                "warning":          res.warning,
                "index_map":        res.index_map,
                "mask":             res.anomaly_mask,
                "preview_b64":      cell_b64,
            })
        except Exception as e:
            app.logger.error(f"{sensor_key} band_depth {mineral} 失败: {e}", exc_info=True)
            sensor_results["band_depth"] = {"error": str(e)}

    # 综合判断:多蚀变叠加(≥2 种异常重叠=高置信)+ 优先级加权评分
    high_pixels = 0
    if masks:
        inter = np.zeros((Hh, Ww), dtype=bool)
        score = np.zeros((Hh, Ww), dtype=np.float32)
        overlap_count = np.zeros((Hh, Ww), dtype=np.int16)
        for t in targets:
            m = masks.get(t["mineral"])
            if m is None:
                continue
            overlap_count += m.astype(np.int16)
            weight = {1: 1.0, 2: 0.5}.get(t["priority"], 0.25)
            score += weight * m.astype(np.float32)
        inter = overlap_count >= 2
        if roi_mask is not None:
            inter &= roi_mask

        if use_structural_prior:
            score = _apply_structural_prior(score, (Hh, Ww), transform, crs, roi_geojson)
            if struct_assoc_out is not None:
                for mineral, a in _structural_assoc_for_sensor(
                        roi_bbox, (Hh, Ww), transform, crs, masks).items():
                    struct_assoc_out.setdefault(mineral, a)

        if basemap_view is not None:
            comp = _render_composite_on_basemap(
                basemap_view, inter, score, transform, crs, roi_geojson, colormap,
                lineaments=struct_lineaments, intersections=struct_intersections,
            )
        else:
            comp = {"intersection": None, "score": None}
        high_pixels = int(inter.sum())
        composites_per_sensor[sensor_key] = {
            "intersection_overlay":   comp.get("intersection"),
            "score_overlay":          comp.get("score"),
            "high_confidence_pixels": high_pixels,
            "shape":                  [Hh, Ww],
            "roi_pixels":             roi_size,
        }
        save_sensors[sensor_key]["composite"] = {
            "intersection_arr":       inter,
            "score_arr":              score,
            "intersection_png":       comp.get("intersection"),
            "score_png":              comp.get("score"),
            "high_confidence_pixels": high_pixels,
        }

    return high_pixels, roi_size


@app.route('/api/analyze_batch', methods=['POST'])
def api_analyze_batch():
    """
    项目驱动的批量蚀变分析。

    body: {
        project_name:    "山东招远庙山金矿4.974km2_1779769677",  必填
        deposit_type:    "斑岩型铜矿",                          必填
        selected_minerals: [...]      (可选,默认全部 priority<=2 且 effective_sensor 有数据的)
        methods:         ["ratio","pca"]                       (默认两种)
        roi_geojson:     {...}        (可选,默认用项目的 .ovkml 全范围)
        threshold_method, k, colormap
    }

    内部:
      - 每个 mineral 自动选 effective_sensor(首选;首选无数据时 fallback)
      - 装载该传感器数据(用 delivery_project 的动态 bn_map)
      - 跑 ratio + pca

    返回 grid: 按矿物 → 传感器 → 方法 三层组织
    """
    data = request.json or {}
    project_name = data.get("project_name", "")
    deposit_type = data.get("deposit_type", "")
    selected_minerals = data.get("selected_minerals")
    methods = data.get("methods") or ["ratio", "pca"]
    roi_geojson = data.get("roi_geojson")
    threshold_method = data.get("threshold_method", "mean_std")
    k = float(data.get("k", 2.0))
    colormap = data.get("colormap", "hot")
    # 可选:用 geo-stru 构造解译产物(距断裂邻近度)对蚀变综合得分做构造控矿加权。
    # 默认关闭,开启后不改变各矿物单独结果,仅重排综合得分(近断裂优先)。
    use_structural_prior = bool(data.get("use_structural_prior", False))
    # 丛林模式:密林区岩石被冠层遮蔽、常规蚀变失效,改为叠加"地植物学胁迫"探测
    # (矿化/微渗漏导致植物受胁迫 → 红边蓝移/叶绿素↓)。默认关闭,不改变现有行为。
    jungle_mode = bool(data.get("jungle_mode", False))

    if not project_name or not deposit_type:
        return jsonify({"error": "缺少 project_name 或 deposit_type"}), 400

    project_dir = DELIVERY_ROOT / project_name
    if not project_dir.is_dir():
        return jsonify({"error": f"项目目录不存在: {project_name}"}), 404

    # 项目可用传感器(目录存在)
    sensors_info = list_available_sensors(project_dir)
    available_sensors = {s["key"] for s in sensors_info}
    if not available_sensors:
        return jsonify({"error": f"项目 {project_name} 冬季子目录无可用传感器数据"}), 404

    # 覆盖校验: 哪些传感器的影像 bbox 与 ROI bbox 有交集(轻量,只读首波段 bounds)
    coverage_info: Dict[str, Dict[str, Any]] = {}
    usable_sensors: set = set()
    roi_bbox = bbox_from_geojson(roi_geojson) if roi_geojson else None

    # 可选:加载 geo-stru 构造解译上下文(断裂线/交汇点/走向统计),供出图叠合与报告。
    # 与 _apply_structural_prior 的评分加权配套;任何缺失都静默降级为 None,不影响分析。
    struct_ctx = None
    if use_structural_prior and roi_bbox:
        try:
            import structural_weighting as _sw
            struct_ctx = _sw.load_structural_context(roi_bbox)
        except Exception as e:
            app.logger.warning(f"加载构造上下文失败,忽略: {e}")
    struct_lineaments = (struct_ctx or {}).get("lineaments")
    struct_intersections = (struct_ctx or {}).get("intersections")
    struct_assoc_by_mineral: Dict[str, Any] = {}
    for sk in available_sensors:
        if roi_bbox is None:
            # 无 ROI 时默认所有可用传感器都参与
            coverage_info[sk] = {"covers": True, "reason": ""}
            usable_sensors.add(sk)
            continue
        ci = check_sensor_coverage(project_dir, sk, roi_bbox)
        coverage_info[sk] = ci
        if ci["covers"]:
            usable_sensors.add(sk)
        else:
            app.logger.info(f"覆盖校验: {sk} 不覆盖 ROI - {ci['reason']}")

    # 多传感器目标视图,自动算 effective_sensor (只在 usable_sensors 中选)
    multi = get_targets_multi_sensor(deposit_type)
    if not multi:
        return jsonify({"error": f"矿床类型 {deposit_type} 无推荐蚀变矿物"}), 400

    # 构建 effective_sensor 表 + 过滤选定矿物
    targets_used = []
    for t in multi:
        preferred = t["preferred_sensor"]
        notes: List[str] = []
        effective = None
        fallback_reason = None

        # 1. 首选传感器是否项目里有数据 + 覆盖 ROI
        if preferred:
            if preferred not in available_sensors:
                fallback_reason = f"首选传感器 {preferred} 在项目中无数据"
            elif preferred not in usable_sensors:
                ci = coverage_info.get(preferred, {})
                fallback_reason = f"首选传感器 {preferred} 影像不覆盖 ROI({ci.get('reason','')})"
            else:
                effective = preferred

        # 2. fallback: 在 usable_sensors 中找首选之外的、且对该矿物有方法的传感器
        if effective is None:
            for s in ("ASTER", "Sentinel2", "Landsat8"):
                if s in usable_sensors:
                    ps = t["per_sensor"].get(s, {})
                    if ps.get("ratio_available") or ps.get("pca_available"):
                        effective = s
                        if fallback_reason:
                            notes.append(f"⚠ {fallback_reason},已 fallback 到 {s}")
                        break

        # 3. 都没找到,记录原因
        if effective is None:
            if fallback_reason:
                # 首选有覆盖问题,fallback 也没有可用方法
                ps_avail = [s for s in ("ASTER","Sentinel2","Landsat8")
                            if s in usable_sensors and
                            (t["per_sensor"].get(s, {}).get("ratio_available") or
                             t["per_sensor"].get(s, {}).get("pca_available"))]
                if not ps_avail:
                    notes.append(f"⚠ {fallback_reason};且其他可用传感器在数据库中均无该矿物的提取方法")
                else:
                    notes.append(f"⚠ {fallback_reason}")
            else:
                notes.append("⚠ 该矿物在项目可用传感器上均无提取方法")

        if selected_minerals is not None and t["mineral"] not in selected_minerals:
            continue
        if selected_minerals is None and t["priority"] > 2:
            continue
        targets_used.append({**t, "effective_sensor": effective, "notes": notes})

    # ─── 丛林模式:注入通用地植物学胁迫探测目标 ───────────────────────
    # 该目标与矿种无关,任何矿床类型在丛林模式下均可叠加。当前走 Sentinel-2 红边波段
    # (B5/B6/B7)的多光谱主路径;EnMAP/PRISMA 仍按各自的 band_depth 正常分析,其红边
    # 胁迫为后续扩展。effective_sensor 选 usable 的 Sentinel-2,缺失则记 note 并置 None。
    if jungle_mode:
        jt_per_sensor = {
            "Sentinel2": {
                "ratio_expr": None, "ratio_available": False,
                "pca_spec": None, "pca_available": False,
                "band_depth_available": False,
                "veg_stress_available": True,
                "veg_stress_spec": {"index": "ndre", "veg_floor": 0.30},
            }
        }
        jungle_target = {
            "mineral":             "丛林地植物胁迫(红边)",
            "zone":                "矿化/微渗漏地表植被冠层",
            "priority":            1,
            "anomaly_type":        "金属毒害/还原胁迫 → 叶绿素↓、红边蓝移",
            "absorption_um":       None,
            "reflectance_peak_um": None,
            "enmap_feature":       None,
            "per_sensor":          jt_per_sensor,
            "preferred_sensor":    "Sentinel2",
            "_methods":            ["veg_stress"],   # 该目标只跑地植物学胁迫法
        }
        jt_eff = "Sentinel2" if "Sentinel2" in usable_sensors else None
        jt_notes = [] if jt_eff else [
            "⚠ 丛林模式地植物学胁迫当前需 Sentinel-2 红边波段,项目无可用 Sentinel-2 覆盖"]
        if selected_minerals is None or jungle_target["mineral"] in (selected_minerals or []):
            targets_used.append({**jungle_target, "effective_sensor": jt_eff, "notes": jt_notes})

    # 按 effective_sensor 分组,每个组装载一次数据
    by_sensor: Dict[str, List[Dict[str, Any]]] = {}
    for t in targets_used:
        if not t["effective_sensor"]:
            continue
        by_sensor.setdefault(t["effective_sensor"], []).append(t)

    # 拉一次卫星底图供所有结果共用(智能优先项目内真彩色,失败 fallback 到 Esri,再失败传感器 RGB)
    basemap_view = None
    if roi_geojson:
        basemap_view = _make_basemap_view(
            roi_geojson, project_dir=project_dir,
            available_sensors=usable_sensors,  # 只用覆盖 ROI 的传感器
        )
        if basemap_view:
            app.logger.info(f"底图就绪 source={basemap_view.get('source')} size={basemap_view['dst_size']}")
        else:
            app.logger.warning("底图全部失败,fallback 到传感器自家真彩色")

    # 跑分析
    per_mineral_results: Dict[str, Dict[str, Any]] = {}
    composites_per_sensor: Dict[str, Dict[str, str]] = {}
    save_sensors: Dict[str, Dict[str, Any]] = {}   # 落盘素材: sensor -> {profile, shape, results, composite}
    total_roi_pixels = 0
    overall_high_pixels = 0

    for sensor_key, tlist in by_sensor.items():
        try:
            image, bn_map, profile = load_sensor_data(project_dir, sensor_key)
        except Exception as e:
            for t in tlist:
                per_mineral_results.setdefault(t["mineral"], {
                    "mineral": t["mineral"], "zone": t["zone"], "priority": t["priority"],
                    "anomaly_type": t["anomaly_type"], "preferred_sensor": t["preferred_sensor"],
                    "effective_sensor": sensor_key, "data_status": "load_failed",
                    "results": {},
                })["results"].setdefault(sensor_key, {"error": f"装载失败: {e}"})
            continue

        if image.max() > 10.0:
            image = image / 10000.0

        H, W = image.shape[1], image.shape[2]
        transform = profile.get("transform") if profile else None
        image_crs = profile.get("crs") if profile else None
        roi_mask = build_roi_mask(
            (H, W), roi_geojson=roi_geojson, transform=transform, image_crs=image_crs
        )

        roi_size = int(roi_mask.sum()) if roi_mask is not None else H * W
        total_roi_pixels = max(total_roi_pixels, roi_size)

        # 收集本传感器的落盘素材(原始 numpy 数组 + 已渲染 PNG,供 save_batch_run 写盘)
        save_sensors[sensor_key] = {
            "profile": profile, "shape": [H, W], "roi_pixels": roi_size,
            "results": [], "composite": None,
        }

        # 每个 mineral 跑两种方法
        sensor_masks: Dict[str, Dict[str, np.ndarray]] = {}  # mineral -> {method: mask}
        vs_index_holder: Dict[str, np.ndarray] = {}          # 捕获地植物学胁迫指数图(丛林模式)

        for t in tlist:
            mineral = t["mineral"]
            entry = per_mineral_results.setdefault(mineral, {
                "mineral":          mineral,
                "zone":             t["zone"],
                "priority":         t["priority"],
                "anomaly_type":     t["anomaly_type"],
                "preferred_sensor": t["preferred_sensor"],
                "effective_sensor": sensor_key,
                "data_status":      "fallback" if (sensor_key != t["preferred_sensor"]) else "ok",
                "notes":            list(t.get("notes") or []),
                "results":          {},
            })
            sensor_results = entry["results"].setdefault(sensor_key, {})
            sensor_masks.setdefault(mineral, {})

            ps = t["per_sensor"].get(sensor_key, {})
            target_spec = {
                "mineral":              mineral,
                "ratio_expr":           ps.get("ratio_expr"),
                "ratio_available":      ps.get("ratio_available", False),
                "pca_spec":             ps.get("pca_spec"),
                "pca_available":        ps.get("pca_available", False),
                "veg_stress_available": ps.get("veg_stress_available", False),
                "veg_stress_spec":      ps.get("veg_stress_spec"),
            }

            # 普通蚀变矿物用全局 methods;丛林地植物胁迫目标用自带 _methods(只 veg_stress)
            t_methods = t.get("_methods") or methods
            for m in t_methods:
                if not ps.get(f"{m}_available"):
                    sensor_results[m] = {"unavailable": True,
                                         "reason": f"{sensor_key} 上 {mineral} 不支持 {m}"}
                    continue
                try:
                    res = analyze_single(
                        image, sensor_key, target_spec, m, roi_mask=roi_mask,
                        threshold_method=threshold_method, k=k, bn_map=bn_map,
                    )
                    if basemap_view is not None:
                        cell_b64 = _render_cell_on_basemap(
                            basemap_view, res.anomaly_mask,
                            profile.get("transform"), profile.get("crs"),
                            roi_geojson, res, sensor_key, colormap,
                        )
                    else:
                        cell_b64 = _render_cell_thumbnail_dyn(image, res, colormap, sensor_key, bn_map)
                    sensor_results[m] = {
                        "anomaly_ratio": round(res.anomaly_ratio * 100, 2),
                        "threshold":     None if not np.isfinite(res.threshold) else float(res.threshold),
                        "warning":       res.warning,
                        "overlay":       f"data:image/png;base64,{cell_b64}" if cell_b64 else None,
                        "ratio_expr":    res.ratio_expr,
                        "pc_used":       res.pc_used,
                        "sign":          res.sign,
                    }
                    sensor_masks[mineral][m] = res.anomaly_mask
                    if m == "veg_stress":
                        vs_index_holder["arr"] = res.index_map
                    save_sensors[sensor_key]["results"].append({
                        "mineral":          mineral,
                        "zone":             t["zone"],
                        "priority":         t["priority"],
                        "anomaly_type":     t["anomaly_type"],
                        "effective_sensor": sensor_key,
                        "data_status":      entry["data_status"],
                        "method":           m,
                        "anomaly_ratio":    round(res.anomaly_ratio * 100, 2),
                        "threshold":        None if not np.isfinite(res.threshold) else float(res.threshold),
                        "ratio_expr":       res.ratio_expr,
                        "pc_used":          res.pc_used,
                        "sign":             res.sign,
                        "warning":          res.warning,
                        "index_map":        res.index_map,
                        "mask":             res.anomaly_mask,
                        "preview_b64":      cell_b64,
                    })
                except Exception as e:
                    app.logger.error(f"analyze {mineral}/{m}/{sensor_key} 失败: {e}", exc_info=True)
                    sensor_results[m] = {"error": str(e)}

        # 该传感器的综合判断
        overall_inter = np.zeros((H, W), dtype=bool)
        score_map = np.zeros((H, W), dtype=np.float32)
        for t in tlist:
            # 丛林地植物胁迫目标不计入"蚀变综合得分":它作为独立证据层(E2)喂给 prospectivity,
            # 计入 score_map 会与 prospectivity 的 veg_stress 层重复计数。
            if t.get("_methods") == ["veg_stress"]:
                continue
            mineral = t["mineral"]
            masks = sensor_masks.get(mineral, {})
            if "ratio" in masks and "pca" in masks:
                inter = masks["ratio"] & masks["pca"]
                overall_inter |= inter
            weight = {1: 1.0, 2: 0.5}.get(t["priority"], 0.25)
            if masks:
                avg = np.mean([m.astype(np.float32) for m in masks.values()], axis=0)
                score_map += weight * avg

        if use_structural_prior:
            score_map = _apply_structural_prior(
                score_map, (H, W), profile.get("transform"), profile.get("crs"), roi_geojson)
            # 异常-构造关联度(逐矿物,基于各方法异常并集)
            mbym = {mineral: _or_masks(mm) for mineral, mm in sensor_masks.items()}
            mbym = {k2: v for k2, v in mbym.items() if v is not None}
            for mineral, a in _structural_assoc_for_sensor(
                    roi_bbox, (H, W), profile.get("transform"), profile.get("crs"), mbym).items():
                struct_assoc_by_mineral.setdefault(mineral, a)

        if basemap_view is not None:
            comp = _render_composite_on_basemap(
                basemap_view, overall_inter, score_map,
                profile.get("transform"), profile.get("crs"),
                roi_geojson, colormap,
                lineaments=struct_lineaments, intersections=struct_intersections,
            )
        else:
            class _FakeBatch: pass
            fb = _FakeBatch()
            fb.sensor = sensor_key
            fb.overall_intersection = overall_inter
            fb.score_map = score_map
            rgb_base = _make_rgb_base_dyn(image, sensor_key, bn_map)
            comp = _render_composite_view_dyn(rgb_base, fb, colormap)
        composites_per_sensor[sensor_key] = {
            "intersection_overlay":   comp.get("intersection"),
            "score_overlay":          comp.get("score"),
            "high_confidence_pixels": int(overall_inter.sum()),
            "shape":                  [H, W],
            "roi_pixels":             roi_size,
        }
        save_sensors[sensor_key]["composite"] = {
            "intersection_arr":       overall_inter,
            "score_arr":              score_map,
            "intersection_png":       comp.get("intersection"),
            "score_png":              comp.get("score"),
            "high_confidence_pixels": int(overall_inter.sum()),
        }
        overall_high_pixels += int(overall_inter.sum())

        # ─── 丛林模式:多证据成矿预测融合(即使蚀变近空也能输出靶区)───
        if jungle_mode:
            try:
                prosp = _compute_prospectivity(
                    score_map, vs_index_holder.get("arr"),
                    (H, W), profile.get("transform"), profile.get("crs"),
                    roi_bbox, roi_mask)
                if prosp is not None:
                    pov = None
                    if basemap_view is not None:
                        rc = _render_composite_on_basemap(
                            basemap_view, prosp.mask, prosp.score,
                            profile.get("transform"), profile.get("crs"),
                            roi_geojson, "inferno",
                            lineaments=struct_lineaments, intersections=struct_intersections)
                        pov = rc.get("score")
                    composites_per_sensor[sensor_key]["prospectivity_overlay"] = pov
                    composites_per_sensor[sensor_key]["prospectivity"] = {
                        "used_layers":        prosp.used_layers,
                        "contributions":      {k: round(v, 3) for k, v in prosp.contributions.items()},
                        "high_target_pixels": int(prosp.mask.sum()),
                        "high_threshold":     None if not np.isfinite(prosp.threshold) else round(float(prosp.threshold), 4),
                    }
                    save_sensors[sensor_key]["composite"]["prospectivity_arr"] = prosp.score
                    save_sensors[sensor_key]["composite"]["prospectivity_png"] = pov
                    app.logger.info(
                        f"丛林模式 prospectivity({sensor_key}): 证据层={prosp.used_layers} "
                        f"靶区像素={int(prosp.mask.sum())}")

                # 天然出露靶向:报告裸露像元占比(出露像元上的常规蚀变才可信)
                _expo_bn = {"Sentinel2": ("B2", "B4", "B8", "B11"),
                            "Landsat8":  ("B2", "B4", "B5", "B6")}.get(sensor_key)
                if _expo_bn and all(b in bn_map for b in _expo_bn):
                    import spectral_unmix as _su
                    bb, rb, nb, sb = (image[bn_map[x]] for x in _expo_bn)
                    expo = _su.detect_natural_exposures(nb, rb, bb, sb, roi_mask=roi_mask)
                    composites_per_sensor[sensor_key]["exposure"] = {
                        "exposure_pixels": expo.stats["exposure_pixels"],
                        "exposure_ratio":  round(expo.stats["exposure_ratio"], 4),
                    }
                    app.logger.info(
                        f"丛林模式 天然出露({sensor_key}): {expo.stats['exposure_pixels']} 像素 "
                        f"({expo.stats['exposure_ratio']*100:.1f}%)")
            except Exception as e:
                app.logger.warning(f"丛林模式 prospectivity 融合失败({sensor_key}): {e}", exc_info=True)

    # ─── 高光谱吸收深度分析(EnMAP / PRISMA,独立于多光谱,作为额外的传感器列)───
    # 两者数据格式不同(EnMAP=单 GeoTIFF+XML,PRISMA=vnir/swir ENVI),但分析路径一致:
    # 装载立方体+波长 → band_depth。装载或分析失败仅告警,不影响其余结果返回。
    for hs_key, hs_loader, hs_aux in (
        (ENMAP_KEY,  load_enmap_data,  get_enmap_metadata_path),
        (PRISMA_KEY, load_prisma_data, get_prisma_metadata_paths),
    ):
        if hs_key not in usable_sensors:
            continue
        if not any(t["per_sensor"].get(hs_key, {}).get("band_depth_available")
                   for t in targets_used):
            continue
        try:
            hs_cube, hs_wl, hs_profile = hs_loader(project_dir)
        except Exception as e:
            app.logger.warning(f"{hs_key} 装载失败,跳过高光谱分析: {e}")
            continue
        aux = hs_aux(project_dir)
        aux_list = aux if isinstance(aux, (list, tuple)) else ([aux] if aux else [])
        high_pixels, hs_roi_size = _run_hyperspectral_analysis(
            hs_key, hs_cube, hs_wl, hs_profile, aux_list,
            targets_used, roi_geojson, basemap_view,
            threshold_method, k, colormap,
            per_mineral_results, composites_per_sensor, save_sensors,
            use_structural_prior=use_structural_prior,
            struct_lineaments=struct_lineaments, struct_intersections=struct_intersections,
            roi_bbox=roi_bbox, struct_assoc_out=struct_assoc_by_mineral,
        )
        total_roi_pixels = max(total_roi_pixels, hs_roi_size)
        overall_high_pixels += high_pixels

    # 也要记录数据缺失/无覆盖的 mineral(effective_sensor=None)
    for t in targets_used:
        if t["effective_sensor"] is None:
            # 状态细分: no_coverage(首选不覆盖且无 fallback)vs no_data(全部传感器无方法)
            notes = t.get("notes") or []
            has_coverage_issue = any("不覆盖" in n or "无数据" in n for n in notes)
            per_mineral_results.setdefault(t["mineral"], {
                "mineral":          t["mineral"],
                "zone":             t["zone"],
                "priority":         t["priority"],
                "anomaly_type":     t["anomaly_type"],
                "preferred_sensor": t["preferred_sensor"],
                "effective_sensor": None,
                "data_status":      "no_coverage" if has_coverage_issue else "no_data",
                "notes":            notes,
                "results":          {},
            })

    # 构造约束汇总(供前端出图说明/报告)+ 逐矿物关联度回填到 grid
    structural_summary = None
    if struct_ctx:
        for mineral, a in struct_assoc_by_mineral.items():
            if mineral in per_mineral_results:
                per_mineral_results[mineral]["structural_assoc"] = a
        structural_summary = {
            "enabled":                   True,
            "aoi_name":                  struct_ctx.get("aoi_name"),
            "n_lineaments":              struct_ctx.get("n_lineaments"),
            "dominant_strikes_deg":      struct_ctx.get("dominant_strikes_deg"),
            "n_intersections":           len(struct_ctx.get("intersections") or []),
            "total_lineament_length_km": struct_ctx.get("total_lineament_length_km"),
            "lineament_density_mean":    struct_ctx.get("lineament_density_mean"),
            "rose_diagram":              _png_to_datauri(struct_ctx.get("rose_diagram_path")),
            "association_by_mineral":    struct_assoc_by_mineral,
        }
    elif use_structural_prior:
        structural_summary = {
            "enabled": True,
            "warning": "未找到与本项目 ROI 相交的构造解译数据,构造约束未生效",
        }

    # 自动落盘(失败只告警,不影响分析返回)
    saved_info = None
    try:
        saved_info = save_batch_run(
            project_name=project_name,
            deposit_type=deposit_type,
            tenant_id=request.headers.get('X-Tenant-Id'),
            params={
                "threshold_method":  threshold_method,
                "k":                 k,
                "colormap":          colormap,
                "methods":           methods,
                "selected_minerals": selected_minerals,
                "use_structural_prior": use_structural_prior,
                "jungle_mode":       jungle_mode,
            },
            roi_geojson=roi_geojson,
            available_sensors=sorted(available_sensors),
            usable_sensors=sorted(usable_sensors),
            coverage_info=coverage_info,
            total_roi_pixels=total_roi_pixels,
            high_confidence_total_pixels=overall_high_pixels,
            sensors=save_sensors,
            structural=structural_summary,
        )
        app.logger.info(f"分析结果已保存: {saved_info['run_id']} "
                        f"({saved_info['n_rasters']} 栅格 / {saved_info['n_previews']} 预览)")
    except Exception as e:
        app.logger.warning(f"保存分析结果失败: {e}", exc_info=True)

    return jsonify({
        "project_name":      project_name,
        "deposit_type":      deposit_type,
        "available_sensors": sorted(available_sensors),
        "usable_sensors":    sorted(usable_sensors),
        "coverage_info":     coverage_info,
        "total_roi_pixels":  total_roi_pixels,
        "grid":              list(per_mineral_results.values()),
        "composites":        composites_per_sensor,
        "high_confidence_total_pixels": overall_high_pixels,
        "structural":        structural_summary,
        "jungle_mode":       jungle_mode,
        "saved":             saved_info,
    })


@app.route('/api/saved_runs', methods=['GET'])
def api_saved_runs():
    """列出已保存的分析 run 摘要(按时间倒序)。?project_name= 可只看某项目。"""
    project = request.args.get("project_name") or None
    try:
        return jsonify({"runs": list_runs(project)})
    except Exception as e:
        app.logger.error(f"列举已保存结果失败: {e}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@app.route('/api/saved_run/<path:run_id>', methods=['GET'])
def api_saved_run(run_id):
    """读取某次 run 的完整 manifest;预览/栅格用相对路径,前端拼 base_url 即可访问。"""
    manifest = load_manifest(run_id)
    if manifest is None:
        return jsonify({"error": "结果不存在"}), 404
    return jsonify({"manifest": manifest, "base_url": "/api/saved_file/" + run_id})


@app.route('/api/saved_file/<path:relpath>', methods=['GET'])
def api_saved_file(relpath):
    """服务 results/ 下的预览 PNG / 下载 GeoTIFF。带越界防护。"""
    target = resolve_file(relpath)
    if target is None:
        return jsonify({"error": "文件不存在或路径非法"}), 404
    return send_file(str(target))


@app.route('/api/download_run_rasters/<path:run_id>', methods=['GET'])
def api_download_run_rasters(run_id):
    """把某次 run 的全部带坐标 GeoTIFF 打成 ZIP 下载(供 GIS 叠加)。带越界防护。"""
    try:
        buf = make_run_rasters_zip(run_id)
    except Exception as e:
        app.logger.error(f"打包栅格 ZIP 失败 {run_id}: {e}", exc_info=True)
        return jsonify({"error": str(e)}), 500
    if buf is None:
        return jsonify({"error": "结果不存在或无可下载的栅格"}), 404
    # 文件名用 run_id 末段(<deposit>_<timestamp>),避免路径分隔符
    safe_tail = run_id.rstrip("/").split("/")[-1] or "rasters"
    return send_file(
        buf, mimetype="application/zip", as_attachment=True,
        download_name=f"{safe_tail}_GeoTIFF_WGS84.zip",
    )


@app.route('/api/analyze_preview', methods=['POST'])
def api_analyze_preview():
    """
    快速预览：只生成指数图缩略图，不做完整统计。
    请求体同 /api/analyze，返回 {"preview": "data:image/png;base64,..."}
    """
    data = request.json or {}
    file_path = data.get("file_path", "")
    sensor_raw = data.get("sensor", "Landsat8/9")
    sensor = sensor_raw.replace("_L2", "")
    mineral = data.get("mineral", "iron_oxide")
    colormap = data.get("colormap", "RdYlGn_r")

    if not file_path or (not os.path.isfile(file_path) and not os.path.isdir(file_path)):
        return jsonify({"error": "文件不存在"}), 400

    try:
        image, _ = read_image(file_path)
        if sensor_raw.endswith("_L2") and image.max() > 10.0:
            image = image / 10000.0

        result = analyze_alteration(image, sensor, mineral)
        b64 = _render_index_map(result.index_map, colormap, title=result.label)
        return jsonify({"preview": f"data:image/png;base64,{b64}" if b64 else None})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ─────────────────────────────────────────────
# 可视化辅助函数
# ─────────────────────────────────────────────

def _render_alteration_overlay(
    image: np.ndarray,
    result,
    colormap: str = "hot",
    dpi: int = 100,
) -> Optional[str]:
    """
    生成蚀变叠加图：底图为 RGB 真彩色，异常区域用 colormap 半透明叠加。
    返回 base64 PNG 字符串，失败返回 None。
    """
    try:
        from alteration_analysis import LANDSAT_BANDS, SENTINEL2_BANDS, ASTER_SWIR_BANDS

        # 底图 RGB
        sensor = result.sensor
        if sensor == "Landsat8/9":
            r_idx, g_idx, b_idx = LANDSAT_BANDS["red"], LANDSAT_BANDS["green"], LANDSAT_BANDS["blue"]
        elif sensor == "Sentinel2":
            r_idx, g_idx, b_idx = SENTINEL2_BANDS["red"], SENTINEL2_BANDS["green"], SENTINEL2_BANDS["blue"]
        else:
            # ASTER SWIR 无真彩色，用前三波段伪彩色
            r_idx, g_idx, b_idx = 2, 1, 0

        n_bands = image.shape[0]
        r_idx = min(r_idx, n_bands - 1)
        g_idx = min(g_idx, n_bands - 1)
        b_idx = min(b_idx, n_bands - 1)

        def _norm(band):
            v = band[np.isfinite(band) & (band > 0)]
            if v.size == 0:
                return np.zeros_like(band)
            lo, hi = np.percentile(v, [2, 98])
            out = (band - lo) / max(hi - lo, 1e-6)
            return np.clip(out, 0, 1)

        rgb = np.stack([_norm(image[r_idx]), _norm(image[g_idx]), _norm(image[b_idx])], axis=2)

        # 绘图
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=dpi)

        # 左：底图 + 异常叠加
        axes[0].imshow(rgb)
        # 掩膜叠加（异常=不透明，正常=透明）
        mask_float = result.anomaly_mask.astype(np.float32)
        cmap = plt.get_cmap(colormap)
        axes[0].imshow(
            np.ma.masked_where(~result.anomaly_mask, mask_float),
            cmap=cmap, alpha=0.6, vmin=0, vmax=1,
            interpolation="nearest"
        )
        axes[0].set_title(f"{result.label} 异常叠加", fontsize=11)
        axes[0].axis("off")

        # 右：连续指数图
        idx = result.index_map.copy()
        idx[~np.isfinite(idx)] = np.nan
        vmin = float(np.nanpercentile(idx, 2))
        vmax = float(np.nanpercentile(idx, 98))
        im = axes[1].imshow(idx, cmap=colormap, vmin=vmin, vmax=vmax, interpolation="nearest")
        axes[1].set_title(f"比值指数图（阈值 {result.threshold:.3f}）", fontsize=11)
        axes[1].axis("off")
        plt.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")

    except Exception:
        return None


def _sensor_rgb_indices(sensor: str, n_bands: int):
    """根据传感器返回 (r, g, b) 在影像中的通道索引,失败回退到前 3 个通道。"""
    from alteration_analysis import LANDSAT_BANDS, SENTINEL2_BANDS
    if sensor == "Landsat8" or sensor == "Landsat8/9":
        idx = (LANDSAT_BANDS["red"], LANDSAT_BANDS["green"], LANDSAT_BANDS["blue"])
    elif sensor == "Sentinel2":
        idx = (SENTINEL2_BANDS["red"], SENTINEL2_BANDS["green"], SENTINEL2_BANDS["blue"])
    else:
        idx = (min(2, n_bands - 1), min(1, n_bands - 1), 0)
    return tuple(min(i, n_bands - 1) for i in idx)


def _make_rgb_base(image: np.ndarray, sensor: str) -> np.ndarray:
    """构建 2-98 百分位拉伸后的 RGB 底图。返回 (H,W,3) float [0,1]。"""
    r_i, g_i, b_i = _sensor_rgb_indices(sensor, image.shape[0])

    def _norm(band):
        v = band[np.isfinite(band) & (band > 0)]
        if v.size == 0:
            return np.zeros_like(band)
        lo, hi = np.percentile(v, [2, 98])
        return np.clip((band - lo) / max(hi - lo, 1e-6), 0, 1)

    return np.stack([_norm(image[r_i]), _norm(image[g_i]), _norm(image[b_i])], axis=2)


def _make_rgb_base_dyn(image: np.ndarray, sensor_key: str, bn_map: Dict[str, int]) -> np.ndarray:
    """
    动态 bn_map 版本: 用波段名直接挑 RGB。
    - ASTER: R=B3N/B3, G=B2, B=B1 (VNIR 真彩近似)
    - Sentinel2: R=B4, G=B3, B=B2
    - Landsat8: R=B4, G=B3, B=B2
    """
    if sensor_key == "ASTER":
        r_name = "B3N" if "B3N" in bn_map else "B3"
        candidates = (r_name, "B2", "B1")
    else:  # Sentinel2 / Landsat8
        candidates = ("B4", "B3", "B2")
    idxs = []
    for n in candidates:
        if n in bn_map:
            idxs.append(bn_map[n])
        else:
            idxs.append(min(len(idxs), image.shape[0] - 1))

    def _norm(band):
        v = band[np.isfinite(band) & (band > 0)]
        if v.size == 0:
            return np.zeros_like(band)
        lo, hi = np.percentile(v, [2, 98])
        return np.clip((band - lo) / max(hi - lo, 1e-6), 0, 1)

    return np.stack([_norm(image[i]) for i in idxs], axis=2)


def _render_cell_thumbnail(
    image: np.ndarray,
    result,
    colormap: str = "hot",
    dpi: int = 90,
    figsize: tuple = (4.2, 4.2),
) -> Optional[str]:
    """
    网格视图的单格缩略图: RGB 底图 + 异常掩膜半透明叠加。
    标题显示方法和异常占比。返回 base64 PNG。
    """
    try:
        rgb = _make_rgb_base(image, result.sensor)
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.imshow(rgb)
        if result.anomaly_mask.any():
            mask_float = result.anomaly_mask.astype(np.float32)
            ax.imshow(
                np.ma.masked_where(~result.anomaly_mask, mask_float),
                cmap=colormap, alpha=0.55, vmin=0, vmax=1, interpolation="nearest",
            )
        if result.method == "pca" and result.pc_used:
            sub = f"PC{result.pc_used} (sign={'+' if result.sign > 0 else '-'})"
        elif result.method == "ratio" and result.ratio_expr:
            sub = result.ratio_expr
        else:
            sub = ""
        title = f"{result.mineral} · {result.method}  {sub}"
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"_render_cell_thumbnail 失败: {e}", exc_info=True)
        return None


def _render_composite_view(image: np.ndarray, batch, colormap: str = "hot") -> Dict[str, Optional[str]]:
    """
    综合判断视图: 返回 {intersection: base64, score: base64}
    intersection: 所有矿物两方法交集的并集
    score: 多蚀变叠加评分(归一化后 colormap 渲染)
    """
    out = {"intersection": None, "score": None}
    rgb = _make_rgb_base(image, batch.sensor)

    # 交集图
    try:
        fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=100)
        ax.imshow(rgb)
        inter = batch.overall_intersection
        if inter is not None and inter.any():
            mask = inter.astype(np.float32)
            ax.imshow(
                np.ma.masked_where(~inter, mask),
                cmap=colormap, alpha=0.7, vmin=0, vmax=1, interpolation="nearest",
            )
        n_pix = int(inter.sum()) if inter is not None else 0
        ax.set_title(f"两方法交集 (高置信度异常 {n_pix} 像素)", fontsize=10)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0); plt.close(fig)
        out["intersection"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite intersection 失败: {e}", exc_info=True)

    # 评分图
    try:
        score = batch.score_map
        if score is None:
            score = np.zeros(image.shape[1:], dtype=np.float32)
        smax = float(score.max())
        score_norm = score / smax if smax > 0 else score
        fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=100)
        ax.imshow(rgb)
        if smax > 0:
            ax.imshow(
                np.ma.masked_where(score_norm < 0.05, score_norm),
                cmap=colormap, alpha=0.65, vmin=0, vmax=1, interpolation="nearest",
            )
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=0, vmax=smax))
            sm.set_array([])
            plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="叠加评分")
        ax.set_title(f"多蚀变叠加评分 (max={smax:.2f})", fontsize=10)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0); plt.close(fig)
        out["score"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite score 失败: {e}", exc_info=True)

    return out


def _render_cell_thumbnail_dyn(
    image: np.ndarray, result, colormap: str, sensor_key: str, bn_map: Dict[str, int],
    dpi: int = 150, figsize: tuple = (7.5, 7.5),
) -> Optional[str]:
    """单格缩略图,但 RGB 底图用动态 bn_map。"""
    try:
        rgb = _make_rgb_base_dyn(image, sensor_key, bn_map)
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.imshow(rgb)
        if result.anomaly_mask.any():
            mask_float = result.anomaly_mask.astype(np.float32)
            ax.imshow(
                np.ma.masked_where(~result.anomaly_mask, mask_float),
                cmap=colormap, alpha=0.55, vmin=0, vmax=1, interpolation="nearest",
            )
        if result.method == "pca" and result.pc_used:
            sub = f"PC{result.pc_used} (sign={'+' if result.sign > 0 else '-'})"
        elif result.method == "ratio" and result.ratio_expr:
            sub = result.ratio_expr
        else:
            sub = ""
        title = f"{result.mineral} · {result.method} · {sensor_key}  {sub}"
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0); plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"_render_cell_thumbnail_dyn 失败: {e}", exc_info=True)
        return None


def _render_composite_view_dyn(rgb_base: np.ndarray, batch_like, colormap: str = "hot") -> Dict[str, Optional[str]]:
    """综合视图,RGB 底图由外部传入(已经 dyn 渲染好)。"""
    out = {"intersection": None, "score": None}
    inter = getattr(batch_like, "overall_intersection", None)
    score = getattr(batch_like, "score_map", None)

    try:
        fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=100)
        ax.imshow(rgb_base)
        if inter is not None and inter.any():
            mask = inter.astype(np.float32)
            ax.imshow(np.ma.masked_where(~inter, mask),
                      cmap=colormap, alpha=0.7, vmin=0, vmax=1, interpolation="nearest")
        n_pix = int(inter.sum()) if inter is not None else 0
        ax.set_title(f"两方法交集 (高置信度 {n_pix} 像素)", fontsize=10)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO(); plt.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0); plt.close(fig)
        out["intersection"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite intersection failed: {e}", exc_info=True)

    try:
        smax = float(score.max()) if score is not None else 0.0
        score_norm = score / smax if (score is not None and smax > 0) else (score if score is not None else None)
        fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=100)
        ax.imshow(rgb_base)
        if score is not None and smax > 0:
            ax.imshow(np.ma.masked_where(score_norm < 0.05, score_norm),
                      cmap=colormap, alpha=0.65, vmin=0, vmax=1, interpolation="nearest")
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=0, vmax=smax))
            sm.set_array([])
            plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="叠加评分")
        ax.set_title(f"多蚀变叠加评分 (max={smax:.2f})", fontsize=10)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO(); plt.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0); plt.close(fig)
        out["score"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite score failed: {e}", exc_info=True)

    return out


# ─────────────────────────────────────────────
# 卫星底图叠加渲染 (v3+)
# ─────────────────────────────────────────────

def _basemap_view_from_local_sensor(
    project_dir: Path,
    sensor_key: str,
    roi_geojson: Dict[str, Any],
    target_max_px: int = 900,
) -> Optional[Dict[str, Any]]:
    """
    用项目内某传感器的 R/G/B 三波段合成真彩色底图,reproject 到 EPSG:4326 ROI bbox 范围。
    Sentinel-2: B4/B3/B2; ASTER: B3N/B2/B1; Landsat8: B4/B3/B2

    若 ROI 在源 raster 上的有效像素数太少(<600 主边),返回 None 让上层降级到 Esri 瓦片
    (高分辨率底图),否则上采样到 target_max_px 即使再大也是糊的。
    """
    try:
        from delivery_project import load_sensor_data
        image, bn_map, profile = load_sensor_data(project_dir, sensor_key)
    except Exception:
        return None

    if sensor_key == "ASTER":
        rgb_names = ("B3N", "B2", "B1")
    else:  # Sentinel2 / Landsat8 都用 B4/B3/B2
        rgb_names = ("B4", "B3", "B2")
    for n in rgb_names:
        if n not in bn_map:
            return None

    src_crs = profile.get("crs")
    src_transform = profile.get("transform")
    if src_crs is None or src_transform is None:
        return None

    # 计算目标:与 ROI bbox 等同,选合适的输出尺寸
    bbox = bbox_from_geojson(roi_geojson)
    if bbox is None:
        return None
    min_lon, min_lat, max_lon, max_lat = bbox
    # 适度向外扩 5% 让 ROI 边界不贴在图边
    pad_lon = (max_lon - min_lon) * 0.05
    pad_lat = (max_lat - min_lat) * 0.05
    pad_bbox = (min_lon - pad_lon, min_lat - pad_lat,
                max_lon + pad_lon, max_lat + pad_lat)

    # 估算源 raster 在 ROI 内的等效像素数:用源像素 GSD 反算
    # src_transform.a 是 X 方向像素大小(度或米,看 CRS),b/d 同理
    try:
        from rasterio.warp import transform_bounds
        # 把 ROI bbox 投到 src_crs,看在源像素坐标下跨多少格
        src_l, src_b, src_r, src_t = transform_bounds("EPSG:4326", src_crs,
                                                      pad_bbox[0], pad_bbox[1],
                                                      pad_bbox[2], pad_bbox[3])
        px_x = abs((src_r - src_l) / src_transform.a)
        px_y = abs((src_t - src_b) / src_transform.e)
        src_max_px = max(px_x, px_y)
        app.logger.info(f"底图 {sensor_key}: ROI 在源 raster 等效像素 {int(px_x)}x{int(px_y)}")
        if src_max_px < 600:
            app.logger.info(f"底图 {sensor_key}: 源像素 {int(src_max_px)} < 600,降级到 Esri")
            return None
    except Exception as e:
        app.logger.warning(f"估算源像素失败: {e}")

    # 输出尺寸:按 ROI 经纬度比例,主边 target_max_px
    lon_extent = pad_bbox[2] - pad_bbox[0]
    lat_extent = pad_bbox[3] - pad_bbox[1]
    aspect = lon_extent / max(lat_extent, 1e-9)
    if aspect >= 1:
        dst_w = target_max_px
        dst_h = max(64, int(target_max_px / aspect))
    else:
        dst_h = target_max_px
        dst_w = max(64, int(target_max_px * aspect))

    from rasterio.transform import from_bounds
    from rasterio.warp import reproject, Resampling
    from rasterio.features import geometry_mask
    dst_transform = from_bounds(pad_bbox[0], pad_bbox[1], pad_bbox[2], pad_bbox[3], dst_w, dst_h)

    bands_out = []
    for n in rgb_names:
        src_band = image[bn_map[n]].astype(np.float32)
        # NaN/无效 → 0,重投影时按 nodata 处理
        src_band = np.where(np.isfinite(src_band), src_band, 0.0)
        dst_band = np.full((dst_h, dst_w), np.nan, dtype=np.float32)
        reproject(
            source=src_band, destination=dst_band,
            src_transform=src_transform, src_crs=src_crs,
            dst_transform=dst_transform, dst_crs="EPSG:4326",
            resampling=Resampling.cubic,
            src_nodata=0.0, dst_nodata=np.nan,
        )
        bands_out.append(dst_band)

    # 2-95 百分位拉伸 + gamma 0.6(提亮暗部,让真彩色卫星图视觉接近常见 RGB 影像)
    def _stretch(b, lo_pct=2, hi_pct=95, gamma=0.6):
        v = b[np.isfinite(b) & (b > 0)]
        if v.size == 0:
            return np.zeros_like(b)
        lo, hi = np.percentile(v, [lo_pct, hi_pct])
        out = np.clip((b - lo) / max(hi - lo, 1e-6), 0, 1)
        out = np.power(out, gamma)
        out[~np.isfinite(b)] = 0
        return out

    rgb = np.stack([_stretch(bands_out[0]), _stretch(bands_out[1]), _stretch(bands_out[2])], axis=2)
    rgb_uint8 = (rgb * 255).astype(np.uint8)

    # ROI 掩膜
    try:
        geom = roi_geojson["geometry"] if "geometry" in roi_geojson else roi_geojson
        roi_mask = geometry_mask([geom], out_shape=(dst_h, dst_w),
                                 transform=dst_transform, invert=True)
    except Exception:
        roi_mask = np.ones((dst_h, dst_w), dtype=bool)

    return {
        "image":     rgb_uint8,
        "bbox":      list(pad_bbox),
        "transform": dst_transform,
        "dst_size":  (dst_h, dst_w),
        "roi_mask":  roi_mask,
        "source":    f"local:{sensor_key}",
    }


def _make_basemap_view(
    roi_geojson: Dict[str, Any],
    target_max_px: int = 1200,
    project_dir: Optional[Path] = None,
    available_sensors: Optional[set] = None,
) -> Optional[Dict[str, Any]]:
    """
    智能选底图:
      1. 项目内 Sentinel-2 真彩色 (优先)
      2. 项目内 ASTER 真彩色
      3. Esri 在线卫星瓦片(兜底)
    """
    # 1+2: 本地传感器优先
    if project_dir and available_sensors:
        for sk in ("Sentinel2", "ASTER", "Landsat8"):
            if sk in available_sensors:
                bm = _basemap_view_from_local_sensor(
                    project_dir, sk, roi_geojson, target_max_px,
                )
                if bm is not None:
                    app.logger.info(f"底图: 使用项目内 {sk} 真彩色")
                    return bm
                else:
                    app.logger.info(f"底图: {sk} 真彩色装载失败,尝试下一个")

    # 3: Esri 兜底
    bbox = bbox_from_geojson(roi_geojson)
    if bbox is None:
        return None
    bm = fetch_satellite_basemap(bbox, target_max_px=target_max_px)
    if bm is None:
        return None
    from rasterio.transform import from_bounds
    from rasterio.features import geometry_mask
    H, W = bm["image"].shape[:2]
    bm_bbox = bm["bbox"]
    transform = from_bounds(bm_bbox[0], bm_bbox[1], bm_bbox[2], bm_bbox[3], W, H)
    try:
        geom = roi_geojson["geometry"] if "geometry" in roi_geojson else roi_geojson
        roi_mask = geometry_mask([geom], out_shape=(H, W), transform=transform, invert=True)
    except Exception:
        roi_mask = np.ones((H, W), dtype=bool)
    app.logger.info(f"底图: fallback 到 Esri 在线 zoom={bm.get('zoom')}")
    return {
        "image":     bm["image"],
        "bbox":      bm_bbox,
        "transform": transform,
        "dst_size":  (H, W),
        "roi_mask":  roi_mask,
        "zoom":      bm.get("zoom"),
        "source":    "esri_online",
    }


def _reproject_mask_to_basemap(
    mask: np.ndarray,
    src_transform,
    src_crs,
    basemap_view: Dict[str, Any],
) -> np.ndarray:
    """把 mask (H,W) bool 从 src CRS reproject 到底图坐标(EPSG:4326)。"""
    from rasterio.warp import reproject, Resampling
    dst_h, dst_w = basemap_view["dst_size"]
    dst = np.zeros((dst_h, dst_w), dtype=np.float32)
    reproject(
        source=mask.astype(np.float32),
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=basemap_view["transform"],
        dst_crs="EPSG:4326",
        resampling=Resampling.nearest,
    )
    return dst > 0.5


def _draw_roi_outline(ax, roi_geojson: Dict[str, Any], basemap_view: Dict[str, Any], color: str = "#10b981"):
    """在 axes 上画 ROI 多边形描边(像素坐标)。"""
    try:
        geom = roi_geojson["geometry"] if "geometry" in roi_geojson else roi_geojson
        coords = geom.get("coordinates")
        if not coords or geom.get("type") not in ("Polygon", "MultiPolygon"):
            return
        rings = coords if geom["type"] == "Polygon" else [r for poly in coords for r in poly]
        from rasterio.transform import rowcol
        for ring in rings:
            xs, ys = [], []
            for lon, lat, *_ in ring:
                row, col = rowcol(basemap_view["transform"], lon, lat)
                xs.append(col); ys.append(row)
            ax.plot(xs, ys, color=color, linewidth=1.8, linestyle="--", alpha=0.9)
    except Exception as e:
        app.logger.warning(f"_draw_roi_outline 失败: {e}")


def _draw_lineaments(ax, lineaments, intersections, basemap_view,
                     line_color: str = "#38bdf8", pt_color: str = "#f59e0b"):
    """在 axes 上叠加地质构造断裂线(实线)与断裂交汇点(圆点)。坐标 EPSG:4326。"""
    try:
        from rasterio.transform import rowcol
        tr = basemap_view["transform"]
        for ln in (lineaments or []):
            xs, ys = [], []
            for lon, lat in ln:
                r, c = rowcol(tr, lon, lat); xs.append(c); ys.append(r)
            if xs:
                ax.plot(xs, ys, color=line_color, linewidth=1.5, alpha=0.95)
        for lon, lat in (intersections or []):
            r, c = rowcol(tr, lon, lat)
            ax.plot([c], [r], marker="o", color=pt_color, markersize=6,
                    markeredgecolor="white", markeredgewidth=0.8, alpha=0.95)
    except Exception as e:
        app.logger.warning(f"_draw_lineaments 失败: {e}")


def _render_cell_on_basemap(
    basemap_view: Dict[str, Any],
    anomaly_mask_src: np.ndarray,    # 保留参数兼容,但不再用于显示
    src_transform,
    src_crs,
    roi_geojson: Dict[str, Any],
    result,
    sensor_key: str,
    colormap: str,
    dpi: int = 150,
    figsize: tuple = (7.5, 7.5),
) -> Optional[str]:
    """
    单格缩略图: 卫星底图 + ROI 内连续蚀变强度的 colormap 渐变叠加。

    把 result.index_map (连续浮点) 用 bilinear 重投影到底图分辨率,
    仅显示 > threshold 的部分,按 (idx - thr) / (max - thr) 归一化到 [0,1]
    用 colormap 渐变着色 — 既平滑又能看出异常强度差异。
    """
    try:
        from rasterio.warp import reproject, Resampling

        idx_src = result.index_map.astype(np.float32)
        # NaN 填 0 以便重采样;后续会用阈值 + ROI 过滤
        idx_for_warp = np.where(np.isfinite(idx_src), idx_src, 0.0)

        dst_h, dst_w = basemap_view["dst_size"]
        idx_dst = np.full((dst_h, dst_w), np.nan, dtype=np.float32)
        reproject(
            source=idx_for_warp,
            destination=idx_dst,
            src_transform=src_transform, src_crs=src_crs,
            dst_transform=basemap_view["transform"], dst_crs="EPSG:4326",
            resampling=Resampling.bilinear,
            src_nodata=0.0, dst_nodata=np.nan,
        )

        thr = float(result.threshold) if np.isfinite(result.threshold) else None
        if thr is None:
            anom_show_any = False
            display = None
        else:
            in_roi = basemap_view["roi_mask"]
            valid = in_roi & np.isfinite(idx_dst)
            # 1. soft 归一化:thr 处 0,ROI 内最大值处 1。低于 thr 的为负值。
            roi_vals = idx_dst[valid]
            if roi_vals.size == 0:
                anom_show_any = False
                display = None
            else:
                vmax = float(np.percentile(roi_vals[roi_vals > thr], 95)) if (roi_vals > thr).any() else float(roi_vals.max())
                denom = max(vmax - thr, 1e-6)
                soft = np.clip((idx_dst - thr) / denom, 0, 1)
                soft = np.where(valid, soft, 0.0)
                # 2. 高斯模糊让边界连续,小噪声被削弱
                from scipy.ndimage import gaussian_filter
                smooth = gaussian_filter(soft, sigma=1.6)
                # 3. ROI 外抹零;低于 0.05 的视为背景透明
                smooth = np.where(in_roi, smooth, 0.0)
                display = np.where(smooth > 0.05, smooth, np.nan)
                anom_show_any = np.isfinite(display).any()

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.imshow(basemap_view["image"])
        if anom_show_any:
            ax.imshow(
                np.ma.masked_invalid(display),
                cmap=colormap, alpha=0.75, vmin=0, vmax=1,
                interpolation="bilinear",
            )
        _draw_roi_outline(ax, roi_geojson, basemap_view)
        if result.method == "pca" and result.pc_used:
            sub = f"PC{result.pc_used} ({'+' if result.sign > 0 else '-'})"
        elif result.method == "ratio" and result.ratio_expr:
            sub = result.ratio_expr
        else:
            sub = ""
        ax.set_title(f"{result.mineral} · {result.method} · {sensor_key}  {sub}", fontsize=9)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0); plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"_render_cell_on_basemap 失败: {e}", exc_info=True)
        return None


def _render_composite_on_basemap(
    basemap_view: Dict[str, Any],
    intersection_src: np.ndarray,
    score_src: np.ndarray,
    src_transform,
    src_crs,
    roi_geojson: Dict[str, Any],
    colormap: str,
    lineaments=None,
    intersections=None,
) -> Dict[str, Optional[str]]:
    """综合判断: 卫星底图 + ROI 内的交集/评分叠加。
    传入 lineaments/intersections 时,叠加地质构造断裂线与交汇点(构造约束综合评分图)。"""
    out = {"intersection": None, "score": None}
    _has_struct = bool(lineaments)
    try:
        inter_dst = _reproject_mask_to_basemap(intersection_src, src_transform, src_crs, basemap_view)
        inter_show = inter_dst & basemap_view["roi_mask"]
        fig, ax = plt.subplots(figsize=(9, 9), dpi=150)
        ax.imshow(basemap_view["image"])
        if inter_show.any():
            # 二值 mask 高斯平滑 → 中心实色,边缘渐淡
            from scipy.ndimage import gaussian_filter
            smooth = gaussian_filter(inter_show.astype(np.float32), sigma=1.8)
            smooth = np.where(basemap_view["roi_mask"], smooth, 0.0)
            # 归一化到 [0, 1](高斯卷积后 max < 1)
            smax = float(smooth.max())
            if smax > 0:
                smooth = smooth / smax
            display = np.where(smooth > 0.08, smooth, np.nan)
            cmap_obj = plt.get_cmap(colormap)
            # 用 colormap 的 0.7~0.95 段,让"高置信度"显眼
            from matplotlib.colors import LinearSegmentedColormap
            colors = [cmap_obj(0.55), cmap_obj(0.95)]
            hot_cmap = LinearSegmentedColormap.from_list("hi_conf", colors)
            ax.imshow(np.ma.masked_invalid(display), cmap=hot_cmap,
                      alpha=0.80, vmin=0, vmax=1, interpolation="bilinear")
        _draw_roi_outline(ax, roi_geojson, basemap_view)
        if _has_struct:
            _draw_lineaments(ax, lineaments, intersections, basemap_view)
        n_pix = int(inter_show.sum())
        ax.set_title(f"两方法交集 (高置信度 {n_pix} 像素)", fontsize=12)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO(); plt.savefig(buf, format="png", bbox_inches="tight", dpi=150)
        buf.seek(0); plt.close(fig)
        out["intersection"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite intersection (basemap) failed: {e}", exc_info=True)

    try:
        # score 是浮点 (H,W),先 reproject 用 bilinear 然后归一化
        from rasterio.warp import reproject, Resampling
        dst_h, dst_w = basemap_view["dst_size"]
        score_dst = np.zeros((dst_h, dst_w), dtype=np.float32)
        reproject(
            source=score_src.astype(np.float32),
            destination=score_dst,
            src_transform=src_transform, src_crs=src_crs,
            dst_transform=basemap_view["transform"], dst_crs="EPSG:4326",
            resampling=Resampling.bilinear,
        )
        # 只在 ROI 内显示
        score_dst[~basemap_view["roi_mask"]] = 0
        smax = float(score_dst.max())
        score_norm = score_dst / smax if smax > 0 else score_dst

        fig, ax = plt.subplots(figsize=(9, 9), dpi=150)
        ax.imshow(basemap_view["image"])
        if smax > 0:
            ax.imshow(np.ma.masked_where(score_norm < 0.05, score_norm),
                      cmap=colormap, alpha=0.65, vmin=0, vmax=1, interpolation="bilinear")
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=0, vmax=smax))
            sm.set_array([])
            plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="叠加评分")
        _draw_roi_outline(ax, roi_geojson, basemap_view)
        if _has_struct:
            _draw_lineaments(ax, lineaments, intersections, basemap_view)
        _title = "构造约束综合评分" if _has_struct else "多蚀变叠加评分"
        ax.set_title(f"{_title} (max={smax:.2f})", fontsize=12)
        ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO(); plt.savefig(buf, format="png", bbox_inches="tight", dpi=150)
        buf.seek(0); plt.close(fig)
        out["score"] = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"composite score (basemap) failed: {e}", exc_info=True)

    return out


def _render_index_map(
    index_map: np.ndarray,
    colormap: str = "RdYlGn_r",
    title: str = "",
    dpi: int = 80,
) -> Optional[str]:
    """生成指数图缩略图，返回 base64 PNG"""
    try:
        idx = index_map.copy().astype(np.float32)
        idx[~np.isfinite(idx)] = np.nan
        fig, ax = plt.subplots(figsize=(4, 4), dpi=dpi)
        vmin = float(np.nanpercentile(idx, 2))
        vmax = float(np.nanpercentile(idx, 98))
        im = ax.imshow(idx, cmap=colormap, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception:
        return None


def _render_overlay_on_bg(
    bg_b64: str,
    result,
    colormap: str = "hot",
    dpi: int = 120,
    alpha: float = 0.6,
) -> Optional[str]:
    """
    将蚀变指数热力图叠加到用户上传的底图上，并在右侧绘制指数图。
    bg_b64: data:image/...;base64,... 格式的字符串
    """
    try:
        from PIL import Image as PILImage

        # 解码底图
        header, encoded = bg_b64.split(",", 1)
        bg_bytes = base64.b64decode(encoded)
        bg_pil = PILImage.open(io.BytesIO(bg_bytes)).convert("RGBA")
        bg_w, bg_h = bg_pil.size

        # 连续指数图归一化到 [0,1]，NaN 和非异常区域设为透明
        idx = result.index_map.copy().astype(np.float32)
        valid = idx[np.isfinite(idx)]
        if valid.size == 0:
            return None
        vmin = float(np.percentile(valid, 2))
        vmax = float(np.percentile(valid, 98))
        idx_norm = np.clip((idx - vmin) / max(vmax - vmin, 1e-6), 0, 1)

        # 只在异常区域显示热力，其余透明
        idx_norm[~result.anomaly_mask] = np.nan

        # 用 colormap 生成 RGBA 热力图
        cmap = plt.get_cmap(colormap)
        heatmap_rgba = cmap(idx_norm)  # (H, W, 4), NaN → (0,0,0,0)
        # NaN 位置 alpha=0，异常位置 alpha=用户指定值
        heatmap_rgba[..., 3] = np.where(np.isfinite(idx_norm), alpha, 0.0)
        heatmap_uint8 = (heatmap_rgba * 255).astype(np.uint8)

        # resize 热力图到底图尺寸
        heat_pil = PILImage.fromarray(heatmap_uint8, mode="RGBA")
        heat_pil = heat_pil.resize((bg_w, bg_h), PILImage.BILINEAR)

        # 合成
        composite = PILImage.alpha_composite(bg_pil, heat_pil).convert("RGB")

        # 双联图：左=底图+热力叠加，右=完整指数图
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=dpi)

        axes[0].imshow(np.array(composite))
        axes[0].set_title(f"{result.label} 热力叠加（底图）", fontsize=11)
        axes[0].axis("off")
        # 添加 colorbar 说明热力含义
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        plt.colorbar(sm, ax=axes[0], fraction=0.046, pad=0.04, label="指数值")

        full_idx = result.index_map.copy().astype(np.float32)
        full_idx[~np.isfinite(full_idx)] = np.nan
        im = axes[1].imshow(full_idx, cmap=colormap, vmin=vmin, vmax=vmax, interpolation="nearest")
        axes[1].set_title(f"完整指数图（阈值 {result.threshold:.3f}）", fontsize=11)
        axes[1].axis("off")
        plt.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")

    except Exception as e:
        app.logger.error(f"_render_overlay_on_bg 失败: {e}", exc_info=True)
        return None


# 用于类型注解（已在顶部导入）


@app.errorhandler(500)
def handle_500(error):
    """错误处理"""
    return jsonify({
        "error": "内部服务器错误",
        "message": str(error)
    }), 500


# ═════════════════════════════════════════════════════════════════════
# Phase 1.5 InSAR 接入 — 与 geo-insar 文件系统订阅式对接
# ═════════════════════════════════════════════════════════════════════

@app.route('/api/insar/stacks', methods=['GET'])
def api_insar_stacks():
    """列出 geo-insar 标准输出目录下所有可分析的 AOI 堆栈。"""
    try:
        import sys as _sys
        if '/opt/deepexplor-services' not in _sys.path:
            _sys.path.insert(0, '/opt/deepexplor-services')
        from commons.insar_broker import scan_available_aois
        stacks = scan_available_aois()
        return jsonify({"stacks": stacks})
    except Exception as e:
        return jsonify({"error": f"扫描 geo-insar 输出失败: {e}"}), 500


@app.route('/api/insar/analyze', methods=['POST'])
def api_insar_analyze():
    """
    针对单个 AOI 跑 InSAR 异常分析。

    Body JSON:
      { "aoi_name": "...", "mode": "coherence_stability"|"velocity_cluster"|"fusion_mineral",
        "velocity_threshold_mm_year": 5.0, "coherence_threshold": 0.3,
        "mineral_anomaly_path": "..."(可选,仅 fusion_mineral 用) }
    """
    try:
        import sys as _sys
        if '/opt/deepexplor-services' not in _sys.path:
            _sys.path.insert(0, '/opt/deepexplor-services')
        from commons.insar_broker import get_stack_path
        from insar_timeseries import load_insar_stack, temporal_velocity_trend
        import insar_analysis as ia

        data = request.get_json() or {}
        aoi = data.get("aoi_name")
        mode = data.get("mode", "coherence_stability")
        if not aoi:
            return jsonify({"error": "缺少 aoi_name"}), 400

        path = get_stack_path(aoi)
        if not path:
            return jsonify({"error": f"找不到 AOI: {aoi}"}), 404

        stack = load_insar_stack(path)
        if stack["n_pairs"] == 0:
            return jsonify({"error": "该 AOI 下没有干涉对"}), 404

        # 先把堆栈拍平成速率图 + 平均相干性
        velocity_info = temporal_velocity_trend(
            stack, min_coherence=float(data.get("coherence_threshold", 0.3))
        )
        velocity = velocity_info["velocity_mm_per_year"]
        if velocity is None:
            return jsonify({"error": "无法计算速率(可能时间基线无效)"}), 500

        # 用所有干涉对的相干性中位数作为代表
        cohs = [c for c in stack["coherences"] if c is not None]
        coherence = np.median(np.stack(cohs, axis=0), axis=0) if cohs else None

        # 调用对应分析
        if mode == "coherence_stability":
            result = ia.coherence_to_stability(
                coherence if coherence is not None else np.ones_like(velocity),
                threshold=float(data.get("coherence_threshold", 0.3)),
            )
        elif mode == "velocity_cluster":
            result = ia.los_velocity_clustering(
                velocity, coherence,
                velocity_threshold_mm_year=float(data.get("velocity_threshold_mm_year", 5.0)),
                coherence_threshold=float(data.get("coherence_threshold", 0.3)),
            )
        elif mode == "fusion_mineral":
            mp = data.get("mineral_anomaly_path")
            if not mp or not os.path.exists(mp):
                return jsonify({"error": "fusion_mineral 模式需要有效的 mineral_anomaly_path"}), 400
            import rasterio as _rio
            with _rio.open(mp) as src:
                mineral = src.read(1).astype(np.float32)
            if mineral.shape != velocity.shape:
                from scipy.ndimage import zoom
                zy = velocity.shape[0] / mineral.shape[0]
                zx = velocity.shape[1] / mineral.shape[1]
                mineral = zoom(mineral, (zy, zx), order=1)
            result = ia.fusion_deformation_mineral(
                velocity, mineral, coherence,
                weight_v=float(data.get("weight_v", 0.5)),
                weight_m=float(data.get("weight_m", 0.5)),
                coherence_threshold=float(data.get("coherence_threshold", 0.3)),
            )
        elif mode == "structural_attribution":
            # 形变构造归因:先聚类活跃形变,再用 geo-stru 距断裂/坡度给每个簇定性
            clusters = ia.los_velocity_clustering(
                velocity, coherence,
                velocity_threshold_mm_year=float(data.get("velocity_threshold_mm_year", 5.0)),
                coherence_threshold=float(data.get("coherence_threshold", 0.3)),
            )
            dist = slp = None
            try:
                import glob as _g, rasterio as _rio, sys as _sys
                if '/opt/deepexplor-services' not in _sys.path:
                    _sys.path.insert(0, '/opt/deepexplor-services')
                import structural_weighting as _sw
                tifs = sorted(_g.glob(os.path.join(path, 'sentinel1_insar', '*', 'los_displacement.tif')))
                bbox = stack["metas"][0].get("aoi_bbox") if stack.get("metas") else None
                if tifs and bbox:
                    with _rio.open(tifs[0]) as s:
                        ref_tr, ref_crs = s.transform, s.crs
                    layers = _sw.load_structural_layers(tuple(bbox), velocity.shape, ref_tr, ref_crs)
                    if layers:
                        dist, slp = layers.get('distance'), layers.get('slope')
            except Exception as _e:
                app.logger.warning(f"加载 geo-stru 构造层失败,归因退化为按速率符号: {_e}")
            result = ia.attribute_deformation(
                clusters.mask, velocity, lineament_distance=dist, slope=slp,
                fault_dist_m=float(data.get("fault_dist_m", 300.0)),
                steep_slope_deg=float(data.get("steep_slope_deg", 15.0)),
            )
        else:
            return jsonify({"error": f"未知 mode: {mode}"}), 400

        # 渲染 PNG(base64)— RdBu_r 双向 cmap 给形变,其他用 result.colormap
        img_b64 = _render_insar_array(
            result.array, result.colormap, result.vmin, result.vmax,
            title=result.name,
        )
        return jsonify({
            "name": result.name,
            "stats": result.stats,
            "unit": result.unit,
            "colormap": result.colormap,
            "vmin": result.vmin,
            "vmax": result.vmax,
            "image_b64": img_b64,
            "velocity_stats": velocity_info["stats"],
        })
    except Exception as e:
        app.logger.error(f"InSAR 分析失败: {e}", exc_info=True)
        return jsonify({"error": str(e), "traceback": traceback.format_exc()}), 500


@app.route('/api/insar/timeseries', methods=['POST'])
def api_insar_timeseries():
    """
    AOI 的时序速率反演 + 相干性衰减建模。
    Body: { "aoi_name": "...", "min_coherence": 0.3 }
    """
    try:
        import sys as _sys
        if '/opt/deepexplor-services' not in _sys.path:
            _sys.path.insert(0, '/opt/deepexplor-services')
        from commons.insar_broker import get_stack_path
        from insar_timeseries import load_insar_stack, temporal_velocity_trend, coherence_decay_model

        data = request.get_json() or {}
        aoi = data.get("aoi_name")
        if not aoi:
            return jsonify({"error": "缺少 aoi_name"}), 400
        path = get_stack_path(aoi)
        if not path:
            return jsonify({"error": f"找不到 AOI: {aoi}"}), 404

        stack = load_insar_stack(path)
        if stack["n_pairs"] == 0:
            return jsonify({"error": "该 AOI 下没有干涉对"}), 404

        velocity_info = temporal_velocity_trend(
            stack, min_coherence=float(data.get("min_coherence", 0.3))
        )
        decay = coherence_decay_model(stack)

        # 速率图渲染
        img_b64 = None
        if velocity_info["velocity_mm_per_year"] is not None:
            img_b64 = _render_insar_array(
                velocity_info["velocity_mm_per_year"],
                cmap_name="RdBu_r",
                title=f"{aoi} LOS 速率 (mm/year)",
                diverging=True,
            )
        return jsonify({
            "aoi_name": aoi,
            "n_pairs": stack["n_pairs"],
            "velocity_stats": velocity_info["stats"],
            "coherence_decay": decay,
            "image_b64": img_b64,
        })
    except Exception as e:
        app.logger.error(f"InSAR 时序分析失败: {e}", exc_info=True)
        return jsonify({"error": str(e), "traceback": traceback.format_exc()}), 500


def _render_insar_array(arr, cmap_name="hot", vmin=None, vmax=None,
                       title="InSAR", diverging=False):
    """渲染 InSAR 数组为 base64 PNG,用 matplotlib(沿用 app.py 现有风格)。"""
    try:
        a = np.asarray(arr, dtype=np.float32)
        if diverging or cmap_name == "RdBu_r":
            finite = a[np.isfinite(a)]
            if finite.size:
                lo = float(np.percentile(finite, 2))
                hi = float(np.percentile(finite, 98))
                mag = max(abs(lo), abs(hi), 1e-6)
                vmin, vmax = -mag, mag
        fig, ax = plt.subplots(figsize=(8, 7))
        fig.patch.set_facecolor('white')
        im = ax.imshow(a, cmap=cmap_name, vmin=vmin, vmax=vmax, origin='lower')
        plt.colorbar(im, ax=ax)
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        buf.seek(0)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        app.logger.error(f"InSAR 渲染失败: {e}", exc_info=True)
        return None


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)
