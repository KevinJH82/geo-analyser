"""
交付数据项目层

约定:
  - 交付根目录: $DELIVERY_ROOT (默认 /Volumes/大硬盘可劲用/DeepExplor/数据留存/交付数据)
  - 项目目录命名: <主名称>_<编号>,如 "山东招远庙山金矿4.974km2_1779769677"
  - 用户上传 ROI 文件,文件名 == 项目目录名,文件后缀通常 .ovkml (奥维) 或 .kml/.geojson
  - 项目目录结构:
      <项目>/data-矿权-冬季（11-3月）/
          ASTER L2/B1.tif B2.tif ... B14.tif
          Sentinel 2 L2/B01.tiff B02.tiff ... B12.tiff B8A.tiff
          (Landsat 8/9 L2 子目录如有则装载)
  - 本次重构只用冬季子目录(植被干扰少)

提供:
  - resolve_project_dir(uploaded_filename) -> Path 或 None
  - parse_roi_file(path) -> GeoJSON Polygon
  - list_available_sensors(project_dir) -> [{key, label, sub_path, n_bands}]
  - load_sensor_data(project_dir, sensor_key) -> (image, bn_to_idx, profile)
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import rasterio


DELIVERY_ROOT = Path(os.environ.get(
    "DELIVERY_ROOT",
    "/Volumes/大硬盘可劲用/DeepExplor/数据留存/交付数据",
))

# 冬季子目录(暂不支持夏季)
WINTER_SUBDIR = "data-矿权-冬季（11-3月）"

# JSON 数据库里用的传感器 key → 项目内子目录"规范名"(仅用于展示/找不到实际目录时的兜底)
SENSOR_DIR_MAP = {
    "ASTER":     "ASTER L2",
    "Sentinel2": "Sentinel 2 L2",
    "Landsat8":  "Landsat 8 L2",
    "Landsat9":  "Landsat 9 L2",
}

# 传感器 key → 目录名前缀(匹配时大小写/分隔符不敏感,且忽略处理级别后缀)。
# 不同交付批次的子目录命名会变,如 "Sentinel 2 L2" / "Sentinel 2 L2A" / "Landsat 8 L2A",
# 这里用前缀匹配避免漏认。Landsat 8 与 Landsat 9 因前缀不同不会互相误配。
SENSOR_DIR_PREFIXES = {
    "ASTER":     ["aster"],
    "Sentinel2": ["sentinel 2", "sentinel2"],
    "Landsat8":  ["landsat 8", "landsat8"],
    "Landsat9":  ["landsat 9", "landsat9"],
    "EnMAP":     ["enmap"],
    "PRISMA":    ["prisma"],
}

# EnMAP(高光谱)是单文件多波段:一个 SPECTRAL_IMAGE.tif(224 波段)+ METADATA.XML,
# 与多光谱的"每波段一文件"不同,需走独立装载路径(load_enmap_data)。
ENMAP_KEY = "EnMAP"
ENMAP_SPECTRAL_FILE = "SPECTRAL_IMAGE.tif"
ENMAP_METADATA_FILE = "METADATA.XML"

# PRISMA(高光谱)以 ENVI 格式交付:vnir / swir 两个二进制立方体 + 同名 .hdr 头文件,
# vnir≈0.40–0.98µm、swir≈0.94–2.50µm,装载时拼成一个立方体(load_prisma_data)。
# 波长、坐标系等信息在 .hdr 里,rasterio(GDAL ENVI 驱动)可直接读取波段 wavelength 标签。
PRISMA_KEY = "PRISMA"
PRISMA_ENVI_FILES = ("vnir", "swir")   # 二进制立方体文件名(同目录有 <name>.hdr)

# 单波段文件名 -> (sort_num, sort_suffix, normalized_name)
# 与 app.py 中 BAND_FILE_PATTERN 兼容,但额外做"去前导 0"归一化
_BAND_RE = re.compile(r"^[Bb](?:and)?(\d+)([A-Za-z]*)\.(tif|tiff)$", re.IGNORECASE)


# ─────────────────────────────────────────────
# 项目定位
# ─────────────────────────────────────────────

def _project_main_name(filename: str) -> str:
    """ovkml/kml/geojson/shp 文件名 → 主名称(去后缀)。"""
    return Path(filename).stem


def resolve_project_dir(uploaded_filename: str) -> Optional[Path]:
    """按上传文件名主名称定位交付目录下的项目目录。"""
    main = _project_main_name(uploaded_filename)
    if not main:
        return None
    candidate = DELIVERY_ROOT / main
    if candidate.is_dir():
        return candidate
    return None


def list_projects() -> List[Dict[str, str]]:
    """列出交付根目录下的所有项目(用于备选下拉)。"""
    if not DELIVERY_ROOT.is_dir():
        return []
    out = []
    for entry in sorted(DELIVERY_ROOT.iterdir()):
        if entry.is_dir() and not entry.name.startswith("."):
            out.append({"name": entry.name, "path": str(entry)})
    return out


# ─────────────────────────────────────────────
# ROI 文件解析
# ─────────────────────────────────────────────

_KML_NS = {"kml": "http://www.opengis.net/kml/2.2"}


def _parse_kml_coordinates(text: str) -> List[List[float]]:
    """从 KML <coordinates> 文本提取 [[lon, lat], ...] 数组。"""
    pts = []
    for token in text.replace("\n", " ").split():
        parts = token.split(",")
        if len(parts) >= 2:
            try:
                pts.append([float(parts[0]), float(parts[1])])
            except ValueError:
                continue
    return pts


def parse_roi_file(path: Path) -> Optional[Dict[str, Any]]:
    """
    解析 ROI 文件,返回 GeoJSON Polygon (EPSG:4326 / CGCS2000 视作等同)。
    支持: .ovkml / .kml (KML 格式)、.geojson / .json
    """
    p = Path(path)
    suffix = p.suffix.lower()

    if suffix in (".ovkml", ".kml"):
        try:
            text = p.read_text(encoding="utf-8")
            # 去 BOM
            if text.startswith("﻿"):
                text = text[1:]
            root = ET.fromstring(text)
            # 兼容命名空间和无命名空间两种情况
            coords_nodes = root.findall(".//kml:Polygon//kml:coordinates", _KML_NS)
            if not coords_nodes:
                coords_nodes = root.findall(".//Polygon//coordinates")
            if not coords_nodes:
                return None
            ring = _parse_kml_coordinates(coords_nodes[0].text or "")
            if len(ring) < 3:
                return None
            # 闭合
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            return {"type": "Polygon", "coordinates": [ring]}
        except Exception:
            return None

    if suffix in (".geojson", ".json"):
        try:
            import json
            data = json.loads(p.read_text(encoding="utf-8"))
            # 兼容 FeatureCollection / Feature / Geometry
            if data.get("type") == "FeatureCollection":
                feats = data.get("features", [])
                if not feats:
                    return None
                return feats[0].get("geometry")
            if data.get("type") == "Feature":
                return data.get("geometry")
            if data.get("type") in ("Polygon", "MultiPolygon"):
                return data
        except Exception:
            return None

    return None


def bbox_from_geojson(geom: Dict[str, Any]) -> Optional[Tuple[float, float, float, float]]:
    """从 GeoJSON Polygon 提取 (minLon, minLat, maxLon, maxLat)。"""
    if not geom:
        return None
    gtype = geom.get("type")
    coords = geom.get("coordinates")
    if not coords:
        return None
    pts = []
    if gtype == "Polygon":
        for ring in coords:
            pts.extend(ring)
    elif gtype == "MultiPolygon":
        for poly in coords:
            for ring in poly:
                pts.extend(ring)
    if not pts:
        return None
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    return min(lons), min(lats), max(lons), max(lats)


# ─────────────────────────────────────────────
# 传感器扫描 + 装载
# ─────────────────────────────────────────────

def _winter_dir(project_dir: Path) -> Optional[Path]:
    """返回项目的冬季子目录;不存在返回 None。"""
    d = project_dir / WINTER_SUBDIR
    return d if d.is_dir() else None


def _norm_dirname(name: str) -> str:
    """目录名归一化:空白/下划线/连字符 → 单空格,小写,便于前缀匹配。"""
    return re.sub(r"[\s_\-]+", " ", name).strip().lower()


def _dir_has_band_files(d: Path) -> bool:
    try:
        return any(_BAND_RE.match(f.name) for f in d.iterdir())
    except OSError:
        return False


def _find_sensor_subdir(winter: Path, sensor_key: str) -> Optional[Path]:
    """
    在冬季目录里按前缀(忽略大小写/分隔符/处理级别后缀如 L2、L2A)找传感器子目录。
    多个候选时优先"含波段文件、名称更短(更接近规范名)"的那个。找不到返回 None。
    """
    prefixes = SENSOR_DIR_PREFIXES.get(sensor_key)
    if not prefixes:
        return None
    candidates = []
    try:
        entries = list(winter.iterdir())
    except OSError:
        return None
    for d in entries:
        if not d.is_dir():
            continue
        norm = _norm_dirname(d.name)
        if any(norm.startswith(p) for p in prefixes):
            candidates.append(d)
    if not candidates:
        return None
    candidates.sort(key=lambda d: (not _dir_has_band_files(d), len(d.name), d.name))
    return candidates[0]


def _find_special_subdir(project_dir: Path, sensor_key: str,
                         is_valid: "callable") -> Optional[Path]:
    """
    定位高光谱传感器(EnMAP/PRISMA)子目录。

    高光谱有时作为单景放在【项目根级】(与季节目录平级,如 "项目/EnMAP L2/"),
    有时又放在季节子目录里(如 "项目/data-...冬季/PRISMA L2/");本函数同时在
    项目根与冬季子目录下查找,并用 is_valid(dir) 判定该目录是否含所需数据文件。
    多候选时取名称更短者。找不到返回 None。
    """
    prefixes = SENSOR_DIR_PREFIXES.get(sensor_key) or []
    if not prefixes:
        return None
    search_roots = [project_dir]
    winter = _winter_dir(project_dir)
    if winter:
        search_roots.append(winter)
    for root in search_roots:
        try:
            entries = list(root.iterdir())
        except OSError:
            continue
        cands = [
            d for d in entries
            if d.is_dir()
            and any(_norm_dirname(d.name).startswith(p) for p in prefixes)
            and is_valid(d)
        ]
        if cands:
            cands.sort(key=lambda d: (len(d.name), d.name))
            return cands[0]
    return None


def _find_enmap_subdir(project_dir: Path) -> Optional[Path]:
    """定位 EnMAP 子目录(需含 SPECTRAL_IMAGE.tif)。"""
    return _find_special_subdir(project_dir, ENMAP_KEY,
                                lambda d: (d / ENMAP_SPECTRAL_FILE).is_file())


def _find_prisma_subdir(project_dir: Path) -> Optional[Path]:
    """定位 PRISMA 子目录(需含 vnir 或 swir ENVI 立方体)。"""
    return _find_special_subdir(
        project_dir, PRISMA_KEY,
        lambda d: any((d / f).is_file() for f in PRISMA_ENVI_FILES),
    )


def list_available_sensors(project_dir: Path) -> List[Dict[str, Any]]:
    """
    扫描项目冬季目录,返回可用传感器列表。
    每项 {key, label, sub_path, n_bands, band_files}
    """
    winter = _winter_dir(project_dir)
    out = []
    for sensor_key in SENSOR_DIR_PREFIXES:
        # 高光谱(EnMAP/PRISMA)可能在项目根级或季节子目录,需独立定位;多光谱在冬季目录内
        if sensor_key == ENMAP_KEY:
            sub = _find_enmap_subdir(project_dir)
        elif sensor_key == PRISMA_KEY:
            sub = _find_prisma_subdir(project_dir)
        elif winter:
            sub = _find_sensor_subdir(winter, sensor_key)
        else:
            sub = None
        if sub is None:
            continue
        # EnMAP:高光谱单文件(SPECTRAL_IMAGE.tif),走独立识别
        if sensor_key == ENMAP_KEY:
            spectral = sub / ENMAP_SPECTRAL_FILE
            if not spectral.is_file():
                continue
            try:
                with rasterio.open(spectral) as src:
                    n_bands = src.count
            except Exception:
                continue
            out.append({
                "key":           sensor_key,
                "label":         sub.name,
                "sub_path":      str(sub),
                "n_bands":       n_bands,
                "band_files":    [ENMAP_SPECTRAL_FILE],
                "hyperspectral": True,
            })
            continue
        # PRISMA:高光谱 ENVI(vnir+swir),波段数为两立方体之和
        if sensor_key == PRISMA_KEY:
            n_bands = 0
            envi_files = []
            for fname in PRISMA_ENVI_FILES:
                f = sub / fname
                if not f.is_file():
                    continue
                try:
                    with rasterio.open(f) as src:
                        n_bands += src.count
                    envi_files.append(fname)
                except Exception:
                    continue
            if not envi_files:
                continue
            out.append({
                "key":           sensor_key,
                "label":         sub.name,
                "sub_path":      str(sub),
                "n_bands":       n_bands,
                "band_files":    envi_files,
                "hyperspectral": True,
            })
            continue
        band_files = sorted(
            [f for f in sub.iterdir() if _BAND_RE.match(f.name)],
            key=lambda f: _sort_key(f.name),
        )
        if not band_files:
            continue
        out.append({
            "key":        sensor_key,
            "label":      sub.name,   # 实际目录名(可能带 L2A 等后缀)
            "sub_path":   str(sub),
            "n_bands":    len(band_files),
            "band_files": [f.name for f in band_files],
        })
    return out


def _sort_key(name: str) -> Tuple[int, str]:
    m = _BAND_RE.match(name)
    if not m:
        return (9999, name)
    return (int(m.group(1)), m.group(2).upper())


def _normalize_bn(name: str) -> Optional[str]:
    """文件名 → JSON 中使用的波段标识。如 B01.tiff → 'B1', B8A.tiff → 'B8A', B3N.tif → 'B3N'。"""
    m = _BAND_RE.match(name)
    if not m:
        return None
    num = int(m.group(1))
    suffix = m.group(2).upper()
    return f"B{num}{suffix}"


def check_sensor_coverage(project_dir: Path, sensor_key: str,
                          roi_bbox_4326: Tuple[float, float, float, float]) -> Dict[str, Any]:
    """
    轻量校验:只读该传感器第一个波段的 bounds,跟 ROI bbox 求交集。
    不装载整景影像,节省时间。

    Returns {
        "covers": bool,
        "sensor_bbox_4326": (minLon, minLat, maxLon, maxLat) | None,
        "reason": str (无覆盖时的解释)
    }
    """
    if sensor_key not in SENSOR_DIR_PREFIXES:
        return {"covers": False, "sensor_bbox_4326": None,
                "reason": f"未知传感器 {sensor_key}"}
    # 高光谱(EnMAP/PRISMA)可能在项目根级,独立定位;多光谱必须在冬季目录内
    if sensor_key == ENMAP_KEY:
        sub = _find_enmap_subdir(project_dir)
    elif sensor_key == PRISMA_KEY:
        sub = _find_prisma_subdir(project_dir)
    else:
        winter = _winter_dir(project_dir)
        if not winter:
            return {"covers": False, "sensor_bbox_4326": None,
                    "reason": "项目无冬季子目录"}
        sub = _find_sensor_subdir(winter, sensor_key)
    if sub is None:
        return {"covers": False, "sensor_bbox_4326": None,
                "reason": f"项目缺少 {sensor_key} 子目录"}

    # 取任意一个波段读 bounds(EnMAP 用 SPECTRAL_IMAGE.tif;PRISMA 用 swir/vnir ENVI)
    if sensor_key == ENMAP_KEY:
        ref_file = sub / ENMAP_SPECTRAL_FILE
        if not ref_file.is_file():
            return {"covers": False, "sensor_bbox_4326": None,
                    "reason": f"{sub.name} 中无 {ENMAP_SPECTRAL_FILE}"}
    elif sensor_key == PRISMA_KEY:
        ref_file = next((sub / f for f in PRISMA_ENVI_FILES if (sub / f).is_file()), None)
        if ref_file is None:
            return {"covers": False, "sensor_bbox_4326": None,
                    "reason": f"{sub.name} 中无 PRISMA ENVI 文件 (vnir/swir)"}
    else:
        band_files = sorted(
            [f for f in sub.iterdir() if _BAND_RE.match(f.name)],
            key=lambda f: _sort_key(f.name),
        )
        if not band_files:
            return {"covers": False, "sensor_bbox_4326": None,
                    "reason": f"{sub.name} 中无波段文件"}
        ref_file = band_files[0]

    try:
        from rasterio.warp import transform_bounds
        with rasterio.open(ref_file) as src:
            sensor_bbox = transform_bounds(src.crs, "EPSG:4326", *src.bounds)
    except Exception as e:
        return {"covers": False, "sensor_bbox_4326": None,
                "reason": f"无法读取影像 bounds: {e}"}

    s_min_lon, s_min_lat, s_max_lon, s_max_lat = sensor_bbox
    r_min_lon, r_min_lat, r_max_lon, r_max_lat = roi_bbox_4326

    overlap = not (s_max_lon < r_min_lon or s_min_lon > r_max_lon
                   or s_max_lat < r_min_lat or s_min_lat > r_max_lat)

    if not overlap:
        reason = (f"{sensor_key} 影像覆盖 [{s_min_lon:.2f}–{s_max_lon:.2f}°E, "
                  f"{s_min_lat:.2f}–{s_max_lat:.2f}°N] 与 ROI "
                  f"[{r_min_lon:.2f}–{r_max_lon:.2f}°E, "
                  f"{r_min_lat:.2f}–{r_max_lat:.2f}°N] 不重叠")
    else:
        reason = ""

    return {"covers": overlap, "sensor_bbox_4326": list(sensor_bbox), "reason": reason}


def load_sensor_data(project_dir: Path, sensor_key: str) -> Tuple[np.ndarray, Dict[str, int], Dict[str, Any]]:
    """
    装载某传感器的所有波段为 (bands, H, W) 数组,同时返回
    BN→通道索引映射(供算法层用)和参考 profile。

    多波段如果分辨率不同,取数量最多的那组,组内裁到最小公共尺寸。
    """
    winter = _winter_dir(project_dir)
    if not winter:
        raise FileNotFoundError(f"项目 {project_dir.name} 缺少冬季子目录")
    if sensor_key not in SENSOR_DIR_PREFIXES:
        raise ValueError(f"不支持的传感器: {sensor_key}")
    sub = _find_sensor_subdir(winter, sensor_key)
    if sub is None:
        raise FileNotFoundError(f"项目 {project_dir.name} 缺少 {sensor_key} 子目录")

    band_files = sorted(
        [f for f in sub.iterdir() if _BAND_RE.match(f.name)],
        key=lambda f: _sort_key(f.name),
    )
    if not band_files:
        raise FileNotFoundError(f"{sub} 中未找到波段文件")

    # 按像元分辨率分组,取主分辨率
    res_groups: Dict[int, List[Path]] = {}
    for bf in band_files:
        with rasterio.open(bf) as src:
            res = round(abs(src.res[0]))
        res_groups.setdefault(res, []).append(bf)
    dominant = max(res_groups, key=lambda r: len(res_groups[r]))
    band_files = res_groups[dominant]

    # 最小公共尺寸
    shapes = {}
    for bf in band_files:
        with rasterio.open(bf) as src:
            shapes[bf] = (src.height, src.width)
    min_h = min(s[0] for s in shapes.values())
    min_w = min(s[1] for s in shapes.values())

    bands = []
    bn_map: Dict[str, int] = {}
    profile = None
    for i, bf in enumerate(band_files):
        bn = _normalize_bn(bf.name)
        if bn is None:
            continue
        with rasterio.open(bf) as src:
            arr = src.read(1)[:min_h, :min_w].astype(np.float32)
            if profile is None:
                profile = src.profile.copy()
        bands.append(arr)
        bn_map[bn] = len(bands) - 1
        # Sentinel-2 容错: 同时注册带前导 0 的别名(如 'B01' → 同 'B1')
        m = _BAND_RE.match(bf.name)
        if m and m.group(2) == "" and bf.name.startswith(("B0", "b0")):
            bn_map[bf.name.split(".")[0].upper()] = bn_map[bn]
        # ASTER 容错: B3N (Nadir 视角) 同时注册为 B3,便于 JSON 表达式直接用 "B3"
        # (后视 B3B 如有,不做 alias,以免覆盖 B3N 注册)
        if m and m.group(2).upper() == "N":
            bare = f"B{int(m.group(1))}"
            bn_map.setdefault(bare, bn_map[bn])

    if not bands:
        raise FileNotFoundError(f"{sub} 中未找到有效波段")

    image = np.stack(bands, axis=0)
    if profile:
        profile.update(count=image.shape[0], height=min_h, width=min_w)
    return image, bn_map, profile or {}


# ─────────────────────────────────────────────
# EnMAP 高光谱装载
# ─────────────────────────────────────────────

def _parse_enmap_metadata(meta_path: Path) -> Dict[str, np.ndarray]:
    """
    解析 EnMAP L2A METADATA.XML,返回 {wavelength_nm, gain, offset, fwhm_nm}(均按波段顺序)。
    缺失时各项可能为 None。
    """
    tree = ET.parse(str(meta_path))
    root = tree.getroot()

    def _strip(tag: str) -> str:
        return tag.split("}")[-1]

    wl: List[float] = []
    gain: List[float] = []
    off: List[float] = []
    fwhm: List[float] = []
    for e in root.iter():
        tag = _strip(e.tag)
        txt = (e.text or "").strip()
        if not txt:
            continue
        try:
            val = float(txt)
        except ValueError:
            continue
        if tag == "wavelengthCenterOfBand":
            wl.append(val)
        elif tag == "GainOfBand":
            gain.append(val)
        elif tag == "OffsetOfBand":
            off.append(val)
        elif tag == "FWHMOfBand":
            fwhm.append(val)

    out: Dict[str, np.ndarray] = {}
    out["wavelength_nm"] = np.array(wl, dtype=np.float64) if wl else None
    out["gain"] = np.array(gain, dtype=np.float64) if gain else None
    out["offset"] = np.array(off, dtype=np.float64) if off else None
    out["fwhm_nm"] = np.array(fwhm, dtype=np.float64) if fwhm else None
    return out


def load_enmap_data(project_dir: Path) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    装载 EnMAP 高光谱数据。

    Returns:
      reflectance: (n_bands, H, W) float32,已乘增益转反射率(0~1),nodata → NaN
      wavelengths_um: (n_bands,) float64,波段中心波长(微米),与 DB 的 absorption_um 同单位
      profile: rasterio profile(含 crs/transform/H/W),供叠加渲染与落盘地理参考用

    EnMAP 是单文件 SPECTRAL_IMAGE.tif(int16, DN)+ METADATA.XML(波长/增益/偏移)。
    """
    sub = _find_enmap_subdir(project_dir)
    if sub is None:
        raise FileNotFoundError(f"项目 {project_dir.name} 缺少 EnMAP 子目录")
    spectral = sub / ENMAP_SPECTRAL_FILE
    if not spectral.is_file():
        raise FileNotFoundError(f"{sub} 中未找到 {ENMAP_SPECTRAL_FILE}")

    with rasterio.open(spectral) as src:
        cube = src.read().astype(np.float32)          # (bands, H, W) DN
        profile = src.profile.copy()
        nodata = src.nodata

    # nodata → NaN
    if nodata is not None:
        cube[cube == nodata] = np.nan

    # 增益/偏移 → 反射率
    meta_path = sub / ENMAP_METADATA_FILE
    wavelengths_um = None
    if meta_path.is_file():
        meta = _parse_enmap_metadata(meta_path)
        wl_nm = meta.get("wavelength_nm")
        gain = meta.get("gain")
        off = meta.get("offset")
        n = cube.shape[0]
        if gain is not None and len(gain) == n:
            cube *= gain.reshape(n, 1, 1).astype(np.float32)
        if off is not None and len(off) == n:
            cube += off.reshape(n, 1, 1).astype(np.float32)
        if wl_nm is not None and len(wl_nm) == n:
            wavelengths_um = (wl_nm / 1000.0)          # nm → μm

    if wavelengths_um is None:
        raise ValueError(
            f"EnMAP 波长信息缺失或与波段数不符({sub / ENMAP_METADATA_FILE}),无法做吸收深度分析"
        )

    profile.update(dtype="float32", count=cube.shape[0], nodata=float("nan"))
    return cube, np.asarray(wavelengths_um, dtype=np.float64), profile


def get_enmap_metadata_path(project_dir: Path) -> Optional[Path]:
    """返回项目 EnMAP 子目录下的 METADATA.XML 路径(存在才返回,否则 None)。

    EnMAP 高光谱的 SPECTRAL_IMAGE.tif 离开 METADATA.XML(波长/增益/偏移)无法解译,
    故下载/落盘 EnMAP 数据时需把该文件一并带上。
    """
    sub = _find_enmap_subdir(project_dir)
    if sub is None:
        return None
    meta_path = sub / ENMAP_METADATA_FILE
    return meta_path if meta_path.is_file() else None


def load_prisma_data(project_dir: Path) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    装载 PRISMA 高光谱数据(ENVI 格式)。

    Returns:
      reflectance:   (n_bands, H, W) float32,L2D 地表反射率(0~1),nodata → NaN
      wavelengths_um:(n_bands,) float64,波段中心波长(微米),与 DB 的 absorption_um 同单位
      profile:       rasterio profile(含 crs/transform/H/W),供叠加渲染与落盘地理参考用

    PRISMA L2 以 vnir(≈0.40–0.98µm)+ swir(≈0.94–2.50µm)两个 ENVI 立方体交付,
    本函数把存在的立方体按"先 vnir 后 swir"拼成一个立方体,波长从各波段的
    GDAL ENVI wavelength 标签读取(单位纳米 → 微米)。
    """
    sub = _find_prisma_subdir(project_dir)
    if sub is None:
        raise FileNotFoundError(f"项目 {project_dir.name} 缺少 PRISMA 子目录")

    parts: List[Tuple[np.ndarray, np.ndarray, Dict[str, Any]]] = []
    for fname in PRISMA_ENVI_FILES:
        f = sub / fname
        if not f.is_file():
            continue
        with rasterio.open(f) as src:
            arr = src.read().astype(np.float32)            # (bands, H, W)
            prof = src.profile.copy()
            nodata = src.nodata
            wl_nm = np.array(
                [float(src.tags(b).get("wavelength", "nan")) for b in range(1, src.count + 1)],
                dtype=np.float64,
            )
        if nodata is not None:
            arr[arr == nodata] = np.nan
        parts.append((arr, wl_nm, prof))

    if not parts:
        raise FileNotFoundError(f"{sub} 中未找到 PRISMA ENVI 文件 (vnir/swir)")

    # 各立方体空间尺寸应一致;保险起见裁到最小公共 H/W 再拼接
    min_h = min(p[0].shape[1] for p in parts)
    min_w = min(p[0].shape[2] for p in parts)
    cube = np.concatenate([p[0][:, :min_h, :min_w] for p in parts], axis=0)
    wl_nm = np.concatenate([p[1] for p in parts], axis=0)

    if not np.all(np.isfinite(wl_nm)):
        raise ValueError(f"PRISMA 波长信息缺失或无法解析({sub}),无法做吸收深度分析")

    wavelengths_um = wl_nm / 1000.0                         # nm → μm
    profile = parts[0][2]                                   # 以 vnir(或首个存在者)为基准
    profile.update(dtype="float32", count=cube.shape[0], height=min_h, width=min_w,
                   nodata=float("nan"))
    return cube, np.asarray(wavelengths_um, dtype=np.float64), profile


def get_prisma_metadata_paths(project_dir: Path) -> List[Path]:
    """返回 PRISMA 的 ENVI 头文件(vnir.hdr / swir.hdr)路径列表。

    PRISMA 的 vnir/swir 二进制立方体离开 .hdr(波长/坐标系/维度)无法解译,
    故下载/落盘 PRISMA 派生结果时把这些头文件一并带上。
    """
    sub = _find_prisma_subdir(project_dir)
    if sub is None:
        return []
    out: List[Path] = []
    for fname in PRISMA_ENVI_FILES:
        hdr = sub / f"{fname}.hdr"
        if hdr.is_file():
            out.append(hdr)
    return out


# ─────────────────────────────────────────────
# 高层入口: 解析上传文件 → 项目元信息
# ─────────────────────────────────────────────

def open_project_from_upload(uploaded_path: Path) -> Optional[Dict[str, Any]]:
    """
    用户上传一个 ROI 文件(任意位置),按主名称定位项目目录,解析 ROI 几何,
    扫描可用传感器,返回完整项目元信息。
    """
    uploaded_path = Path(uploaded_path)
    project_dir = resolve_project_dir(uploaded_path.name)
    if not project_dir:
        return {"error": f"在交付目录中未找到与 '{_project_main_name(uploaded_path.name)}' 同名的项目",
                "delivery_root": str(DELIVERY_ROOT)}

    roi_geom = parse_roi_file(uploaded_path)
    if not roi_geom:
        return {"error": f"无法解析 ROI 文件 {uploaded_path.name}(支持 .ovkml/.kml/.geojson)",
                "project_dir": str(project_dir)}

    sensors = list_available_sensors(project_dir)
    bbox = bbox_from_geojson(roi_geom)

    return {
        "project_name":    project_dir.name,
        "project_dir":     str(project_dir),
        "roi_geojson":     roi_geom,
        "bbox":            list(bbox) if bbox else None,
        "available_sensors": sensors,
        "winter_subdir":   WINTER_SUBDIR,
    }


if __name__ == "__main__":
    print(f"交付根: {DELIVERY_ROOT}")
    print(f"可用项目: {[p['name'] for p in list_projects()]}")
    print()

    test_kml = (DELIVERY_ROOT / "山东招远庙山金矿4.974km2_1779769677"
                / "山东招远庙山金矿4.974km2_1779769677.ovkml")
    if test_kml.exists():
        info = open_project_from_upload(test_kml)
        print(f"项目: {info['project_name']}")
        print(f"BBox: {info['bbox']}")
        print(f"ROI 顶点数: {len(info['roi_geojson']['coordinates'][0])}")
        print("可用传感器:")
        for s in info["available_sensors"]:
            print(f"  - {s['key']:<10} ({s['label']}) {s['n_bands']} 波段: {s['band_files']}")

        # 测装载
        if info["available_sensors"]:
            sk = info["available_sensors"][0]["key"]
            img, bn_map, prof = load_sensor_data(Path(info["project_dir"]), sk)
            print(f"\n{sk} 装载: {img.shape}, CRS={prof.get('crs')}")
            print(f"BN→idx: {bn_map}")
