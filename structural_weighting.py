"""
structural_weighting.py — 用 geo-stru 构造解译产物增强蚀变异常分析

提供两类可选增强(默认不启用,不改变现有分析行为):

1. 地形/光照归一化(terrain_normalize):用 geo-stru 的山体阴影压制
   "向阳坡反射率高、背阳坡低"导致的坡度假异常(对应 geo-analyser
   "异常检测不考虑地形"短板)。

2. 构造控矿加权(apply_structural_weighting):用"距断裂距离"把分散的
   蚀变像元约束为沿构造展布的连贯靶区——成矿受构造控制,断裂(尤其交汇点)
   附近的蚀变才是真正有利区。

产物发现复用 commons/structural_broker(bbox 相交 + metadata.json),
栅格按目标分析格网重投影对齐。全部为只读消费,不修改 geo-stru 输出。
"""

import sys
import numpy as np
from typing import Dict, Optional, Tuple

_REPO = "/opt/deepexplor-services"
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


def _reproject_to(path: str, ref_shape, ref_transform, ref_crs) -> Optional[np.ndarray]:
    """把一个 GeoTIFF 重投影/重采样到目标分析格网。失败返回 None。"""
    try:
        import rasterio
        from rasterio.warp import reproject, Resampling
        dst = np.full(ref_shape, np.nan, dtype=np.float32)
        with rasterio.open(path) as src:
            reproject(
                source=rasterio.band(src, 1), destination=dst,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=ref_transform, dst_crs=ref_crs,
                resampling=Resampling.bilinear,
            )
        return dst
    except Exception:
        return None


def load_structural_layers(
    bbox: Tuple[float, float, float, float],
    ref_shape, ref_transform, ref_crs,
) -> Optional[Dict[str, np.ndarray]]:
    """
    发现与 bbox 相交的 geo-stru 产物,把所需图层对齐到分析格网。

    Returns
    -------
    {'distance': arr|None, 'density': arr|None, 'hillshade': arr|None,
     'slope': arr|None, 'aoi_name': str} 或 None(无产物)
    """
    try:
        from commons.structural_broker import find_structural_for_bbox, get_product_path
    except Exception:
        return None

    matches = find_structural_for_bbox(bbox)
    if not matches:
        return None
    entry = matches[0]  # 取第一个相交 AOI(同区一般唯一)

    out = {'aoi_name': entry.get('aoi_name', '')}
    for key, pkey in [('distance', 'distance_to_lineament'),
                      ('density', 'lineament_density'),
                      ('hillshade', 'hillshade'),
                      ('slope', 'slope')]:
        path = get_product_path(entry, pkey)
        out[key] = _reproject_to(path, ref_shape, ref_transform, ref_crs) if path else None
    return out


def proximity_from_distance(distance: np.ndarray, scale_m: Optional[float] = None) -> np.ndarray:
    """距断裂距离(米)→ 邻近度[0,1](指数衰减,近断裂=1)。"""
    valid = np.isfinite(distance)
    if not valid.any():
        return np.zeros_like(distance, dtype=np.float32)
    if scale_m is None:
        scale_m = float(np.nanmedian(distance[valid])) or 1.0
    prox = np.where(valid, np.exp(-distance / (scale_m + 1e-9)), 0.0)
    return prox.astype(np.float32)


def terrain_normalize(index_map: np.ndarray, hillshade: np.ndarray) -> np.ndarray:
    """
    用相对光照(山体阴影)归一化蚀变指数,压制坡向/光照引起的亮度差异。
    hillshade 任意尺度;内部归一到均值=1 后做除法,再保持原量纲。
    """
    hs = hillshade.astype(np.float32)
    finite = np.isfinite(hs) & (hs > 0)
    if not finite.any():
        return index_map
    rel = hs / (np.nanmean(hs[finite]) + 1e-6)
    rel = np.clip(rel, 0.3, 3.0)  # 限幅,避免极端阴影区放大噪声
    out = index_map.astype(np.float32) / rel
    out[~np.isfinite(out)] = np.nan
    return out


def apply_structural_weighting(
    anomaly: np.ndarray,
    proximity: np.ndarray,
    weight: float = 0.15,
) -> np.ndarray:
    """
    构造控矿加权:near-fault 上调、far-fault 轻度下调。
    weight=0 时还原为原异常;默认 0.15。
    """
    p = np.clip(np.nan_to_num(proximity, nan=0.0), 0, 1)
    return (anomaly * (1.0 - weight + weight * p)).astype(np.float32)


# ─────────────────────────────────────────────
# 断裂矢量 / 交汇点 / 异常-构造关联度(供出图叠合与报告标注)
# 复用 commons/structural_broker 的 bbox 相交发现;只读消费,不改 geo-stru 输出。
# ─────────────────────────────────────────────

def _segment_intersections(lines, max_points: int = 300):
    """纯 numpy 两两线段求交,返回交点 [(lon,lat),...];成矿最有利的断裂交汇部位。"""
    segs = []
    for ln in lines:
        for i in range(len(ln) - 1):
            segs.append((ln[i], ln[i + 1]))
    pts = []
    n = len(segs)
    for i in range(n):
        (x1, y1), (x2, y2) = segs[i]
        for j in range(i + 1, n):
            (x3, y3), (x4, y4) = segs[j]
            d = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3)
            if abs(d) < 1e-12:
                continue  # 平行/共线
            t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / d
            u = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / d
            if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
                px, py = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
                if all(abs(px - qx) > 1e-9 or abs(py - qy) > 1e-9 for qx, qy in pts):
                    pts.append((float(px), float(py)))
                    if len(pts) >= max_points:
                        return pts
    return pts


def load_structural_context(bbox) -> Optional[Dict]:
    """
    发现与 bbox 相交的 geo-stru 构造产物,返回断裂线矢量/交汇点/走向统计/玫瑰图路径。
    供出图叠合与报告使用。无匹配/无矢量返回 None。
    """
    try:
        from commons.structural_broker import find_structural_for_bbox, get_product_path
    except Exception:
        return None
    matches = find_structural_for_bbox(bbox) if bbox else []
    if not matches:
        return None
    entry = matches[0]
    stats = entry.get("structural_stats", {}) or {}

    lines = []
    gj_path = get_product_path(entry, "lineaments_geojson")
    if gj_path:
        try:
            import json
            with open(gj_path, "r", encoding="utf-8") as f:
                gj = json.load(f)
            for feat in gj.get("features", []):
                geom = (feat or {}).get("geometry") or {}
                gtype = geom.get("type")
                if gtype == "LineString":
                    lines.append([(float(x), float(y)) for x, y in geom["coordinates"]])
                elif gtype == "MultiLineString":
                    for part in geom["coordinates"]:
                        lines.append([(float(x), float(y)) for x, y in part])
        except Exception:
            pass

    return {
        "aoi_name":                  entry.get("aoi_name", ""),
        "n_lineaments":              int(stats.get("n_lineaments", len(lines))),
        "dominant_strikes_deg":      stats.get("dominant_strikes_deg", []) or [],
        "total_lineament_length_km": stats.get("total_lineament_length_km"),
        "lineament_density_mean":    stats.get("lineament_density_mean"),
        "lineaments":                lines,                       # [[(lon,lat),...], ...]
        "intersections":            _segment_intersections(lines),  # [(lon,lat), ...]
        "rose_diagram_path":        get_product_path(entry, "rose_diagram"),
    }


def lineament_association(anomaly_mask, distance, buffer_m: float = 300.0) -> Optional[Dict]:
    """
    异常-构造关联度:异常像元到最近断裂的中位距离(米)+ 落在 buffer 内的比例。
    distance 为对齐到分析格网的'距断裂距离'数组(load_structural_layers 提供);
    distance/掩膜缺失返回 None。
    """
    if distance is None or anomaly_mask is None:
        return None
    m = np.asarray(anomaly_mask, dtype=bool)
    n = int(m.sum())
    if n == 0:
        return {"n_anomaly": 0, "median_dist_m": None, "frac_within_buffer": None}
    d = np.asarray(distance, dtype=np.float32)
    vals = d[m & np.isfinite(d)]
    if vals.size == 0:
        return {"n_anomaly": n, "median_dist_m": None, "frac_within_buffer": None}
    return {
        "n_anomaly":          n,
        "median_dist_m":      float(np.median(vals)),
        "frac_within_buffer": float(np.mean(vals <= buffer_m)),
        "buffer_m":           float(buffer_m),
    }
