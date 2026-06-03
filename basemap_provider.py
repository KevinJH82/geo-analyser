"""
在线卫星瓦片底图提供器

按 ROI 经纬度 bbox 自动选 zoom,从 Esri World Imagery (公共,免 token) 拉瓦片拼合,
裁到 bbox,返回 (RGB ndarray, actual_bbox_lonlat, transform)。

公共瓦片 URL 模板:
  Esri:        https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}
  OSM(备用):  https://tile.openstreetmap.org/{z}/{x}/{y}.png

约定:
  - bbox 用 (min_lon, min_lat, max_lon, max_lat),EPSG:4326
  - 输出 RGB 是 (H, W, 3) uint8;transform 是仿射变换(WebMercator EPSG:3857)
  - 网络/瓦片失败时返回 None,由调用方 fallback
"""

from __future__ import annotations

import io
import math
import os
import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_TILE_SOURCE = os.environ.get(
    "BASEMAP_TILE_URL",
    "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
)
DEFAULT_USER_AGENT = "deepexplor-geo-analyser/1.0"

# 单次请求超时(秒)与每瓦片重试
TILE_TIMEOUT = 8
TILE_RETRIES = 2
# 最大瓦片数(避免 ROI 过大撑爆下载)
MAX_TILES = 64
TILE_PX = 256


# ─────────────────────────────────────────────
# WebMercator 坐标 ↔ 瓦片 xyz
# ─────────────────────────────────────────────

def _lonlat_to_tile(lon: float, lat: float, z: int) -> Tuple[float, float]:
    """返回浮点瓦片坐标 (x, y);整数部分=瓦片号,小数部分=瓦片内位置。"""
    lat_rad = math.radians(lat)
    n = 2.0 ** z
    x = (lon + 180.0) / 360.0 * n
    y = (1.0 - math.log(math.tan(lat_rad) + 1.0 / math.cos(lat_rad)) / math.pi) / 2.0 * n
    return x, y


def _tile_to_lonlat(x: float, y: float, z: int) -> Tuple[float, float]:
    n = 2.0 ** z
    lon = x / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * y / n)))
    return lon, math.degrees(lat_rad)


def _choose_zoom(bbox: Tuple[float, float, float, float], target_tiles: int = 16) -> int:
    """选一个 zoom 让 bbox 内瓦片数接近 target_tiles(默认 4×4=16)。"""
    min_lon, min_lat, max_lon, max_lat = bbox
    for z in range(2, 20):
        x0, y0 = _lonlat_to_tile(min_lon, max_lat, z)
        x1, y1 = _lonlat_to_tile(max_lon, min_lat, z)
        n_tiles = max(1, math.ceil(x1) - math.floor(x0)) * max(1, math.ceil(y1) - math.floor(y0))
        if n_tiles >= target_tiles:
            return z
    return 18


# ─────────────────────────────────────────────
# 瓦片下载与拼合
# ─────────────────────────────────────────────

def _fetch_tile(url_template: str, z: int, x: int, y: int):
    try:
        import requests
    except ImportError:
        return None
    url = url_template.format(z=z, x=x, y=y)
    for attempt in range(TILE_RETRIES + 1):
        try:
            r = requests.get(url, timeout=TILE_TIMEOUT,
                             headers={"User-Agent": DEFAULT_USER_AGENT})
            r.raise_for_status()
            from PIL import Image as PILImage
            return PILImage.open(io.BytesIO(r.content)).convert("RGB")
        except Exception as e:
            if attempt == TILE_RETRIES:
                logger.warning(f"瓦片 {z}/{x}/{y} 下载失败: {e}")
                return None


def fetch_satellite_basemap(
    bbox: Tuple[float, float, float, float],
    target_max_px: int = 1024,
    tile_url: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    主入口: 按 ROI bbox(WGS84 lonlat) 拉卫星瓦片合成。

    Returns dict {
        "image":          (H, W, 3) uint8 RGB,
        "bbox":           实际覆盖 bbox (min_lon, min_lat, max_lon, max_lat),
        "transform":      rasterio 仿射变换(WebMercator EPSG:3857),
        "crs":            "EPSG:3857",
        "zoom":           int,
    } 或 None(失败)。
    """
    if tile_url is None:
        tile_url = DEFAULT_TILE_SOURCE

    min_lon, min_lat, max_lon, max_lat = bbox
    if max_lon <= min_lon or max_lat <= min_lat:
        return None

    z = _choose_zoom(bbox, target_tiles=16)

    # 瓦片 xy 范围 (注意 y 是从北到南递增)
    x0_f, y0_f = _lonlat_to_tile(min_lon, max_lat, z)
    x1_f, y1_f = _lonlat_to_tile(max_lon, min_lat, z)
    x0, x1 = math.floor(x0_f), math.ceil(x1_f)
    y0, y1 = math.floor(y0_f), math.ceil(y1_f)
    n_tiles = (x1 - x0) * (y1 - y0)
    if n_tiles > MAX_TILES:
        # 降 zoom 让瓦片数变少
        while n_tiles > MAX_TILES and z > 2:
            z -= 1
            x0_f, y0_f = _lonlat_to_tile(min_lon, max_lat, z)
            x1_f, y1_f = _lonlat_to_tile(max_lon, min_lat, z)
            x0, x1 = math.floor(x0_f), math.ceil(x1_f)
            y0, y1 = math.floor(y0_f), math.ceil(y1_f)
            n_tiles = (x1 - x0) * (y1 - y0)

    # 下载并拼接
    from PIL import Image as PILImage
    big = PILImage.new("RGB", ((x1 - x0) * TILE_PX, (y1 - y0) * TILE_PX), (40, 40, 40))
    any_ok = False
    for xi in range(x0, x1):
        for yi in range(y0, y1):
            tile = _fetch_tile(tile_url, z, xi, yi)
            if tile is not None:
                big.paste(tile, ((xi - x0) * TILE_PX, (yi - y0) * TILE_PX))
                any_ok = True
    if not any_ok:
        return None

    # 实际拼合后的 bbox(对齐到瓦片边界,通常比 ROI bbox 略大)
    actual_min_lon, actual_max_lat = _tile_to_lonlat(x0, y0, z)
    actual_max_lon, actual_min_lat = _tile_to_lonlat(x1, y1, z)

    # 裁切回 ROI bbox
    px_per_lon = big.size[0] / (actual_max_lon - actual_min_lon)
    px_per_lat = big.size[1] / (actual_max_lat - actual_min_lat)
    left   = int(round((min_lon - actual_min_lon) * px_per_lon))
    right  = int(round((max_lon - actual_min_lon) * px_per_lon))
    top    = int(round((actual_max_lat - max_lat) * px_per_lat))
    bottom = int(round((actual_max_lat - min_lat) * px_per_lat))
    left, top = max(0, left), max(0, top)
    right, bottom = min(big.size[0], right), min(big.size[1], bottom)
    if right > left and bottom > top:
        big = big.crop((left, top, right, bottom))
        cropped_bbox = (min_lon, min_lat, max_lon, max_lat)
    else:
        cropped_bbox = (actual_min_lon, actual_min_lat, actual_max_lon, actual_max_lat)

    # 控制最大边
    H, W = big.size[1], big.size[0]
    scale = min(target_max_px / max(H, W), 1.0)
    if scale < 1.0:
        big = big.resize((int(W * scale), int(H * scale)), PILImage.BILINEAR)

    return {
        "image":  np.array(big),
        "bbox":   cropped_bbox,
        "zoom":   z,
        "source": tile_url,
    }


if __name__ == "__main__":
    # demo: 招远的 bbox
    bbox = (120.42938167, 37.17746472, 120.45604833, 37.20024250)
    print("拉招远 bbox 卫星底图...")
    res = fetch_satellite_basemap(bbox)
    if res is None:
        print("失败(网络?)")
    else:
        print(f"图像 shape: {res['image'].shape}, zoom={res['zoom']}")
        print(f"bbox: {res['bbox']}")
        from PIL import Image as PILImage
        PILImage.fromarray(res["image"]).save("/tmp/basemap_test.png")
        print("已保存到 /tmp/basemap_test.png")
