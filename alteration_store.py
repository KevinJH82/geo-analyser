"""
蚀变分析结果本地持久化层。

每次 /api/analyze_batch 成功后自动落盘到应用本地 results/ 目录,用于:
  1. 后续查看 / 历史对比;
  2. 作为基础数据供其他子系统(clustering_analysis / insar_analysis 等)消费。

目录结构:
    results/
      <project_safe>/
        <deposit_safe>_<YYYYMMDD-HHMMSS>/        一次运行 = 一个 run 目录
          manifest.json                          唯一权威索引(见 save_batch_run)
          rasters/
            <sensor>__<mineral>__<method>__index.tif   float32, nodata=NaN, 带 crs+transform
            <sensor>__<mineral>__<method>__mask.tif     uint8 0/1
            <sensor>__composite__intersection.tif       uint8
            <sensor>__composite__score.tif              float32
          previews/
            <sensor>__<mineral>__<method>.png           复用已渲染的 base64,不重渲染
            <sensor>__composite__intersection.png
            <sensor>__composite__score.png
        latest.json                              {deposit_type: 最新 run_id},供子系统免扫描取最新

交接格式说明(给其他子系统):
    rasters/ 下的 GeoTIFF 即基础数据。*__index.tif 是连续蚀变指数(float32),
    *__mask.tif 是阈值掩膜(uint8 0/1),均带 crs+transform,rasterio.open 直接读,
    正是 clustering_analysis / insar_analysis 期望的 index_map / anomaly_mask 数组。
"""

from __future__ import annotations

import io
import os
import re
import json
import base64
import shutil
import zipfile
import logging
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

# 结果根目录:默认应用本地 results/,允许 RESULTS_ROOT 环境变量覆盖(部署时可指向共享盘)
RESULTS_ROOT = Path(os.environ.get(
    "RESULTS_ROOT",
    str(Path(__file__).parent / "results"),
))

# 同一 (项目+矿床类型) 自动保存会累积,保留最近 N 个 run,删更旧的
MAX_RUNS_PER_DEPOSIT = int(os.environ.get("MAX_RUNS_PER_DEPOSIT", "30"))

SCHEMA_VERSION = "1.0"

_SAFE_RE = re.compile(r"[^\w一-龥-]+")


def results_root() -> Path:
    return RESULTS_ROOT


def _safe_name(s: Any) -> str:
    """与前端 reportFilename 同款清洗:非[字/数/汉字/-]→_,去首尾下划线。"""
    return _SAFE_RE.sub("_", str(s if s is not None else "")).strip("_") or "x"


# ─────────────────────────────────────────────
# 写盘原语
# ─────────────────────────────────────────────

def _write_raster(path: Path, array: np.ndarray, profile: Optional[Dict[str, Any]],
                  dtype: str, nodata: Optional[float]) -> None:
    """写单波段 GeoTIFF。profile 带 crs/transform 时附地理参考,否则写普通 GeoTIFF。"""
    import rasterio
    H, W = array.shape
    prof: Dict[str, Any] = {
        "driver": "GTiff", "height": H, "width": W, "count": 1,
        "dtype": dtype, "compress": "deflate",
    }
    if nodata is not None:
        prof["nodata"] = nodata
    if profile:
        if profile.get("crs") is not None:
            prof["crs"] = profile["crs"]
        if profile.get("transform") is not None:
            prof["transform"] = profile["transform"]
    with rasterio.open(str(path), "w", **prof) as dst:
        dst.write(array.astype(dtype), 1)


def _write_png_b64(path: Path, data: Optional[str]) -> bool:
    """把 base64 PNG(可带 data:image/png;base64, 前缀)写为文件。返回是否写成功。"""
    if not data:
        return False
    if data.startswith("data:"):
        data = data.split(",", 1)[1]
    try:
        path.write_bytes(base64.b64decode(data))
        return True
    except Exception as e:
        logger.warning("写 PNG 预览失败 %s: %s", path, e)
        return False


def _crs_to_str(crs: Any) -> Optional[str]:
    if crs is None:
        return None
    try:
        return crs.to_string()
    except Exception:
        return str(crs)


def _transform_to_list(transform: Any) -> Optional[List[float]]:
    if transform is None:
        return None
    try:
        return [float(x) for x in list(transform)[:6]]
    except Exception:
        return None


# ─────────────────────────────────────────────
# 主入口:保存一次批量分析
# ─────────────────────────────────────────────

def save_batch_run(
    project_name: str,
    deposit_type: str,
    params: Dict[str, Any],
    roi_geojson: Optional[Dict[str, Any]],
    available_sensors: List[str],
    usable_sensors: List[str],
    coverage_info: Dict[str, Any],
    total_roi_pixels: int,
    high_confidence_total_pixels: int,
    sensors: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    把一次批量分析落盘并返回 {run_id, run_dir, n_rasters, n_previews}。

    sensors: {
      <sensor_key>: {
        "profile": <rasterio profile, 含 crs/transform>,
        "shape":   [H, W],
        "roi_pixels": int,
        "results": [ {mineral, zone, priority, anomaly_type, effective_sensor,
                      data_status, method, anomaly_ratio, threshold, ratio_expr,
                      pc_used, sign, warning,
                      index_map: np.ndarray, mask: np.ndarray, preview_b64: str|None} ],
        "composite": {intersection_arr: np, score_arr: np,
                      intersection_png: str|None, score_png: str|None,
                      high_confidence_pixels: int} | None,
      }
    }
    """
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d-%H%M%S")
    project_safe = _safe_name(project_name)
    deposit_safe = _safe_name(deposit_type)

    run_id = f"{project_safe}/{deposit_safe}_{timestamp}"   # 相对 RESULTS_ROOT 的标识
    run_dir = RESULTS_ROOT / run_id
    raster_dir = run_dir / "rasters"
    preview_dir = run_dir / "previews"
    raster_dir.mkdir(parents=True, exist_ok=True)
    preview_dir.mkdir(parents=True, exist_ok=True)

    n_rasters = 0
    n_previews = 0
    n_aux = 0
    manifest_sensors: Dict[str, Any] = {}
    manifest_results: List[Dict[str, Any]] = []
    manifest_composites: Dict[str, Any] = {}

    for sensor_key, sdata in sensors.items():
        profile = sdata.get("profile") or {}
        manifest_sensors[sensor_key] = {
            "crs":        _crs_to_str(profile.get("crs")),
            "transform":  _transform_to_list(profile.get("transform")),
            "shape":      sdata.get("shape"),
            "roi_pixels": sdata.get("roi_pixels"),
        }

        # 随源数据携带的辅助文件(如 EnMAP 的 METADATA.XML — 离开它高光谱无法解译)
        # 复制到 rasters/<sensor>__<原名>,与该传感器的栅格放在一起,一并进 ZIP
        aux_rels: List[str] = []
        for aux in (sdata.get("aux_files") or []):
            src = Path(aux)
            if not src.is_file():
                logger.warning("辅助文件不存在,跳过: %s", src)
                continue
            rel = f"rasters/{sensor_key}__{src.name}"
            try:
                shutil.copyfile(src, run_dir / rel)
                aux_rels.append(rel)
                n_aux += 1
            except Exception as e:
                logger.warning("复制辅助文件失败 %s: %s", src, e)
        if aux_rels:
            manifest_sensors[sensor_key]["aux_files"] = aux_rels

        # 逐个 (矿物 × 方法) 结果
        for r in sdata.get("results", []):
            base = f"{sensor_key}__{_safe_name(r['mineral'])}__{r['method']}"
            entry: Dict[str, Any] = {
                "mineral":          r.get("mineral"),
                "zone":             r.get("zone"),
                "priority":         r.get("priority"),
                "anomaly_type":     r.get("anomaly_type"),
                "effective_sensor": r.get("effective_sensor"),
                "data_status":      r.get("data_status"),
                "sensor":           sensor_key,
                "method":           r.get("method"),
                "anomaly_ratio":    r.get("anomaly_ratio"),
                "threshold":        r.get("threshold"),
                "ratio_expr":       r.get("ratio_expr"),
                "pc_used":          r.get("pc_used"),
                "sign":             r.get("sign"),
                "warning":          r.get("warning"),
                "index_tif":        None,
                "mask_tif":         None,
                "preview_png":      None,
            }
            try:
                idx = r.get("index_map")
                if idx is not None:
                    rel = f"rasters/{base}__index.tif"
                    _write_raster(run_dir / rel, np.asarray(idx, dtype=np.float32),
                                  profile, "float32", float("nan"))
                    entry["index_tif"] = rel
                    n_rasters += 1
                msk = r.get("mask")
                if msk is not None:
                    rel = f"rasters/{base}__mask.tif"
                    _write_raster(run_dir / rel, np.asarray(msk).astype(np.uint8),
                                  profile, "uint8", None)
                    entry["mask_tif"] = rel
                    n_rasters += 1
                rel_png = f"previews/{base}.png"
                if _write_png_b64(run_dir / rel_png, r.get("preview_b64")):
                    entry["preview_png"] = rel_png
                    n_previews += 1
            except Exception as e:
                logger.warning("保存结果栅格失败 %s: %s", base, e)
            manifest_results.append(entry)

        # 综合判断
        comp = sdata.get("composite")
        if comp:
            cmeta: Dict[str, Any] = {
                "high_confidence_pixels": comp.get("high_confidence_pixels"),
                "intersection_tif": None, "score_tif": None,
                "intersection_png": None, "score_png": None,
            }
            try:
                inter = comp.get("intersection_arr")
                if inter is not None:
                    rel = f"rasters/{sensor_key}__composite__intersection.tif"
                    _write_raster(run_dir / rel, np.asarray(inter).astype(np.uint8),
                                  profile, "uint8", None)
                    cmeta["intersection_tif"] = rel
                    n_rasters += 1
                score = comp.get("score_arr")
                if score is not None:
                    rel = f"rasters/{sensor_key}__composite__score.tif"
                    _write_raster(run_dir / rel, np.asarray(score, dtype=np.float32),
                                  profile, "float32", None)
                    cmeta["score_tif"] = rel
                    n_rasters += 1
                rel_png = f"previews/{sensor_key}__composite__intersection.png"
                if _write_png_b64(run_dir / rel_png, comp.get("intersection_png")):
                    cmeta["intersection_png"] = rel_png
                    n_previews += 1
                rel_png = f"previews/{sensor_key}__composite__score.png"
                if _write_png_b64(run_dir / rel_png, comp.get("score_png")):
                    cmeta["score_png"] = rel_png
                    n_previews += 1
            except Exception as e:
                logger.warning("保存综合栅格失败 %s: %s", sensor_key, e)
            manifest_composites[sensor_key] = cmeta

    manifest = {
        "schema_version":               SCHEMA_VERSION,
        "run_id":                       run_id,
        "created_at":                   now.isoformat(timespec="seconds"),
        "project_name":                 project_name,
        "deposit_type":                 deposit_type,
        "params":                       params,
        "roi_geojson":                  roi_geojson,
        "available_sensors":            available_sensors,
        "usable_sensors":               usable_sensors,
        "coverage_info":                coverage_info,
        "total_roi_pixels":             total_roi_pixels,
        "high_confidence_total_pixels": high_confidence_total_pixels,
        "sensors":                      manifest_sensors,
        "results":                      manifest_results,
        "composites":                   manifest_composites,
    }
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 更新 deposit 级 latest 指针,供子系统免扫描取最新
    _update_latest(project_safe, deposit_type, run_id)

    # 保留策略:超出上限删最旧 run
    try:
        prune_old_runs(project_safe, deposit_safe, keep=MAX_RUNS_PER_DEPOSIT)
    except Exception as e:
        logger.warning("清理旧 run 失败: %s", e)

    return {
        "run_id":      run_id,
        "run_dir":     str(run_dir),
        "n_rasters":   n_rasters,
        "n_previews":  n_previews,
        "n_aux":       n_aux,
    }


def _update_latest(project_safe: str, deposit_type: str, run_id: str) -> None:
    latest_path = RESULTS_ROOT / project_safe / "latest.json"
    data: Dict[str, str] = {}
    if latest_path.is_file():
        try:
            data = json.loads(latest_path.read_text(encoding="utf-8"))
        except Exception:
            data = {}
    data[deposit_type] = run_id
    latest_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def prune_old_runs(project_safe: str, deposit_safe: str, keep: int = 30) -> int:
    """同一 (项目+矿床) 只保留最近 keep 个 run,返回删除数量。"""
    if keep <= 0:
        return 0
    proj_dir = RESULTS_ROOT / project_safe
    if not proj_dir.is_dir():
        return 0
    prefix = f"{deposit_safe}_"
    # run 目录名以时间戳结尾,字典序即时间序
    runs = sorted(
        [d for d in proj_dir.iterdir() if d.is_dir() and d.name.startswith(prefix)],
        key=lambda d: d.name,
    )
    to_delete = runs[:-keep] if len(runs) > keep else []
    for d in to_delete:
        shutil.rmtree(d, ignore_errors=True)
    if to_delete:
        logger.info("prune_old_runs: %s/%s 删除 %d 个旧 run(保留 %d)",
                    project_safe, deposit_safe, len(to_delete), keep)
    return len(to_delete)


# ─────────────────────────────────────────────
# 读取 / 列举
# ─────────────────────────────────────────────

def _manifest_summary(manifest: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "run_id":         manifest.get("run_id"),
        "project_name":   manifest.get("project_name"),
        "deposit_type":   manifest.get("deposit_type"),
        "created_at":     manifest.get("created_at"),
        "n_results":      len(manifest.get("results", [])),
        "usable_sensors": manifest.get("usable_sensors", []),
        "high_confidence_total_pixels": manifest.get("high_confidence_total_pixels"),
    }


def list_runs(project_name: Optional[str] = None) -> List[Dict[str, Any]]:
    """列出已保存的 run 摘要,按时间倒序。project_name 给定时仅该项目。"""
    if not RESULTS_ROOT.is_dir():
        return []
    if project_name:
        proj_dirs = [RESULTS_ROOT / _safe_name(project_name)]
    else:
        proj_dirs = [d for d in RESULTS_ROOT.iterdir() if d.is_dir()]

    out: List[Dict[str, Any]] = []
    for pdir in proj_dirs:
        if not pdir.is_dir():
            continue
        for run in pdir.iterdir():
            mf = run / "manifest.json"
            if not mf.is_file():
                continue
            try:
                manifest = json.loads(mf.read_text(encoding="utf-8"))
                out.append(_manifest_summary(manifest))
            except Exception as e:
                logger.warning("读取 manifest 失败 %s: %s", mf, e)
    out.sort(key=lambda r: r.get("created_at") or "", reverse=True)
    return out


def load_manifest(run_id: str) -> Optional[Dict[str, Any]]:
    """按 run_id(相对路径)读取完整 manifest;不存在返回 None。"""
    # 规范化并做越界防护
    root = RESULTS_ROOT.resolve()
    mf = (RESULTS_ROOT / run_id / "manifest.json").resolve()
    if not (mf == root or str(mf).startswith(str(root) + os.sep)):
        return None
    if not mf.is_file():
        return None
    try:
        return json.loads(mf.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning("读取 manifest 失败 %s: %s", mf, e)
        return None


def resolve_file(relpath: str) -> Optional[Path]:
    """把相对路径解析为 RESULTS_ROOT 内的真实文件,越界或不存在返回 None。"""
    root = RESULTS_ROOT.resolve()
    target = (RESULTS_ROOT / relpath).resolve()
    if not (target == root or str(target).startswith(str(root) + os.sep)):
        return None
    return target if target.is_file() else None


def resolve_run_dir(run_id: str) -> Optional[Path]:
    """把 run_id(相对路径)解析为 RESULTS_ROOT 内的真实 run 目录,越界或不存在返回 None。"""
    root = RESULTS_ROOT.resolve()
    target = (RESULTS_ROOT / run_id).resolve()
    if not (target == root or str(target).startswith(str(root) + os.sep)):
        return None
    return target if target.is_dir() else None


def _reproject_tif_bytes(path: Path, dst_epsg: int = 4326) -> bytes:
    """
    把一个 GeoTIFF 重投影到目标坐标系(默认 EPSG:4326 经纬度),返回新 GeoTIFF 字节。
    源无 crs 或已是目标 crs 时,原样返回文件字节(不重采样)。
    掩膜/uint8 用最近邻,连续值(float)用双线性,保留 nodata。
    """
    import rasterio
    from rasterio.crs import CRS
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    dst_crs = CRS.from_epsg(dst_epsg)
    with rasterio.open(path) as src:
        if src.crs is None or src.crs == dst_crs:
            return path.read_bytes()

        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds,
        )
        profile = src.profile.copy()
        profile.update(
            driver="GTiff", crs=dst_crs, transform=transform,
            width=width, height=height, compress="deflate",
        )
        resampling = (Resampling.nearest if src.dtypes[0].startswith("uint")
                      else Resampling.bilinear)

        with rasterio.io.MemoryFile() as mem:
            with mem.open(**profile) as dst:
                for b in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, b),
                        destination=rasterio.band(dst, b),
                        src_transform=src.transform, src_crs=src.crs,
                        dst_transform=transform, dst_crs=dst_crs,
                        src_nodata=src.nodata, dst_nodata=src.nodata,
                        resampling=resampling,
                    )
            return mem.read()


def make_run_rasters_zip(run_id: str, dst_epsg: Optional[int] = 4326) -> Optional[io.BytesIO]:
    """
    把某次 run 的带地理参考 GeoTIFF(rasters/ 下全部 .tif)打成内存 ZIP,
    供前端"下载全部影像"一次性下载,可直接在 QGIS/ArcGIS 中作为栅格图层叠加。

    dst_epsg 给定时(默认 4326 经纬度 WGS84)把每个栅格重投影到该坐标系,
    使任何软件都能读到经纬度位置信息;传 None 则保留原始(常为 UTM 投影)坐标系。

    ZIP 内结构:
        rasters/<sensor>__<mineral>__<method>__index.tif   连续蚀变指数(float32)
        rasters/<sensor>__<mineral>__<method>__mask.tif    阈值掩膜(uint8 0/1)
        rasters/<sensor>__composite__intersection.tif       综合交集(uint8)
        rasters/<sensor>__composite__score.tif               综合评分(float32)
        rasters/EnMAP__METADATA.XML                          EnMAP 波长/增益元数据(随源携带)
        manifest.json                                        索引(含 crs/transform/参数)

    所有 .tif 均带 crs+transform;辅助文件(如 EnMAP METADATA.XML)原样打包不做投影。
    无可打包栅格时返回 None。
    """
    run_dir = resolve_run_dir(run_id)
    if run_dir is None:
        return None
    raster_dir = run_dir / "rasters"
    tifs = sorted(raster_dir.glob("*.tif")) if raster_dir.is_dir() else []
    if not tifs:
        return None

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for tif in tifs:
            if dst_epsg is not None:
                try:
                    zf.writestr(f"rasters/{tif.name}", _reproject_tif_bytes(tif, dst_epsg))
                    continue
                except Exception as e:
                    logger.warning("重投影 %s 失败(改用原始坐标系): %s", tif.name, e)
            zf.write(tif, arcname=f"rasters/{tif.name}")
        # 辅助文件(如 EnMAP METADATA.XML):非栅格,原样打包不做投影
        for aux in sorted(raster_dir.glob("*")):
            if not aux.is_file() or aux.suffix.lower() in (".tif", ".tiff"):
                continue
            zf.write(aux, arcname=f"rasters/{aux.name}")
        manifest = run_dir / "manifest.json"
        if manifest.is_file():
            zf.write(manifest, arcname="manifest.json")
    buf.seek(0)
    return buf
