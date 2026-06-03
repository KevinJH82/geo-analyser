"""
干扰剔除模块 - 通过多种指数计算生成掩膜，去除植被、水体、建筑物、云等非地质信息干扰
支持 Landsat 8/9、Sentinel-2 等多光谱遥感数据
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
import warnings

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# 波段配置
# ─────────────────────────────────────────────

@dataclass
class BandConfig:
    """传感器波段索引配置（0-based）"""
    blue:  int
    green: int
    red:   int
    nir:   int
    swir1: int
    swir2: int
    name:  str = "custom"


LANDSAT8_BANDS = BandConfig(
    blue=1, green=2, red=3, nir=4, swir1=5, swir2=6, name="Landsat8/9"
)

SENTINEL2_BANDS = BandConfig(
    blue=1, green=2, red=3, nir=7, swir1=9, swir2=10, name="Sentinel-2"
    # 读入顺序: B01(0),B02(1),B03(2),B04(3),B05(4),B06(5),B07(6),B08(7),B8A(8),B11(9),B12(10)
)


# ─────────────────────────────────────────────
# 指数计算
# ─────────────────────────────────────────────

def _safe_divide(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """避免除零，返回 float32"""
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(b != 0, a.astype(np.float32) / b.astype(np.float32), 0.0)
    return result.astype(np.float32)


def calc_ndvi(nir: np.ndarray, red: np.ndarray) -> np.ndarray:
    """
    归一化植被指数 NDVI = (NIR - Red) / (NIR + Red)
    范围 [-1, 1]，植被通常 > 0.2
    """
    return _safe_divide(nir - red, nir + red)


def calc_ndwi(green: np.ndarray, nir: np.ndarray) -> np.ndarray:
    """
    归一化水体指数 NDWI = (Green - NIR) / (Green + NIR)
    范围 [-1, 1]，水体通常 > 0.0
    """
    return _safe_divide(green - nir, green + nir)


def calc_mndwi(green: np.ndarray, swir1: np.ndarray) -> np.ndarray:
    """
    改进归一化水体指数 MNDWI = (Green - SWIR1) / (Green + SWIR1)
    对混浊水体效果优于 NDWI
    """
    return _safe_divide(green - swir1, green + swir1)


def calc_ndbi(swir1: np.ndarray, nir: np.ndarray) -> np.ndarray:
    """
    归一化建筑指数 NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)
    建筑/裸地通常 > 0.0
    """
    return _safe_divide(swir1 - nir, swir1 + nir)


def calc_bsi(blue: np.ndarray, red: np.ndarray,
             nir: np.ndarray, swir1: np.ndarray) -> np.ndarray:
    """
    裸土指数 BSI = ((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))
    用于区分裸露岩石/土壤与植被
    """
    numerator   = (swir1 + red) - (nir + blue)
    denominator = (swir1 + red) + (nir + blue)
    return _safe_divide(numerator, denominator)


# ─────────────────────────────────────────────
# 掩膜生成
# ─────────────────────────────────────────────

@dataclass
class MaskThresholds:
    """各掩膜阈值（可按研究区调整）"""
    ndvi_veg:    float = 0.20   # NDVI > 此值 → 植被
    ndwi_water:  float = 0.00   # NDWI > 此值 → 水体（初判）
    mndwi_water: float = 0.00   # MNDWI > 此值 → 水体（确认）
    ndbi_build:  float = 0.10   # NDBI > 此值 且 NDVI < 0.1 → 建筑
    cloud_blue:  float = 0.25   # 蓝波段反射率 > 此值 → 疑似云（需归一化到[0,1]）
    snow_ndsi:   float = 0.40   # NDSI > 此值 → 雪/冰


def make_vegetation_mask(ndvi: np.ndarray,
                         threshold: float = 0.20) -> np.ndarray:
    """植被掩膜：True = 植被（需剔除）"""
    return ndvi > threshold


def make_water_mask(ndwi: np.ndarray, mndwi: np.ndarray,
                    ndwi_thr: float = 0.00,
                    mndwi_thr: float = 0.00) -> np.ndarray:
    """水体掩膜：True = 水体（需剔除），NDWI 与 MNDWI 联合判断减少误判"""
    return (ndwi > ndwi_thr) | (mndwi > mndwi_thr)


def make_buildup_mask(ndbi: np.ndarray, ndvi: np.ndarray,
                      ndbi_thr: float = 0.10) -> np.ndarray:
    """建筑物掩膜：True = 建筑/城镇（需剔除）"""
    return (ndbi > ndbi_thr) & (ndvi < 0.10)


def make_cloud_mask(blue: np.ndarray, nir: np.ndarray,
                    blue_thr: float = 0.25) -> np.ndarray:
    """
    简单云掩膜（适用于已归一化到 [0,1] 的反射率数据）：
    高蓝波段反射率 + 高 NIR → 厚云
    """
    return (blue > blue_thr) & (nir > blue_thr * 0.8)


def make_snow_mask(green: np.ndarray, swir1: np.ndarray,
                   ndsi_thr: float = 0.40) -> np.ndarray:
    """
    雪/冰掩膜 NDSI = (Green - SWIR1) / (Green + SWIR1)
    """
    ndsi = _safe_divide(green - swir1, green + swir1)
    return ndsi > ndsi_thr


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

@dataclass
class RemovalResult:
    """干扰剔除结果"""
    combined_mask:    np.ndarray          # True = 受干扰像元（剔除）
    vegetation_mask:  np.ndarray
    water_mask:       np.ndarray
    buildup_mask:     np.ndarray
    cloud_mask:       np.ndarray
    snow_mask:        np.ndarray
    indices: dict     = field(default_factory=dict)   # 各指数数组
    stats:   dict     = field(default_factory=dict)   # 统计信息


def remove_interference(
    image: np.ndarray,
    band_cfg: BandConfig = LANDSAT8_BANDS,
    thresholds: Optional[MaskThresholds] = None,
    nodata: float = 0.0,
) -> RemovalResult:
    """
    对多光谱影像执行完整的干扰剔除流程。

    Parameters
    ----------
    image : np.ndarray
        shape = (bands, rows, cols)，反射率值需归一化到 [0, 1]
    band_cfg : BandConfig
        波段索引配置
    thresholds : MaskThresholds
        各指数判断阈值
    nodata : float
        无效值标识

    Returns
    -------
    RemovalResult
    """
    if thresholds is None:
        thresholds = MaskThresholds()

    rows, cols = image.shape[1], image.shape[2]

    # 提取各波段（float32）
    def band(idx: int) -> np.ndarray:
        return image[idx].astype(np.float32)

    blue  = band(band_cfg.blue)
    green = band(band_cfg.green)
    red   = band(band_cfg.red)
    nir   = band(band_cfg.nir)
    swir1 = band(band_cfg.swir1)

    # 无效像元掩膜（所有波段均为 nodata 视为无效）
    nodata_mask = np.all(image == nodata, axis=0)

    # ── 计算指数 ──────────────────────────────
    ndvi  = calc_ndvi(nir, red)
    ndwi  = calc_ndwi(green, nir)
    mndwi = calc_mndwi(green, swir1)
    ndbi  = calc_ndbi(swir1, nir)
    bsi   = calc_bsi(blue, red, nir, swir1)

    # ── 生成各类掩膜 ──────────────────────────
    veg_mask   = make_vegetation_mask(ndvi, thresholds.ndvi_veg)
    water_mask = make_water_mask(ndwi, mndwi,
                                  thresholds.ndwi_water,
                                  thresholds.mndwi_water)
    build_mask = make_buildup_mask(ndbi, ndvi, thresholds.ndbi_build)
    cloud_mask = make_cloud_mask(blue, nir, thresholds.cloud_blue)
    snow_mask  = make_snow_mask(green, swir1, thresholds.snow_ndsi)

    # ── 合并掩膜 ─────────────────────────────
    combined = (veg_mask | water_mask | build_mask |
                cloud_mask | snow_mask | nodata_mask)

    # ── 统计 ─────────────────────────────────
    total = rows * cols
    def pct(m): return float(np.sum(m)) / total * 100

    stats = {
        "total_pixels":      total,
        "vegetation_pct":    round(pct(veg_mask),   2),
        "water_pct":         round(pct(water_mask),  2),
        "buildup_pct":       round(pct(build_mask),  2),
        "cloud_pct":         round(pct(cloud_mask),  2),
        "snow_pct":          round(pct(snow_mask),   2),
        "combined_pct":      round(pct(combined),    2),
        "valid_geology_pct": round(100 - pct(combined), 2),
    }

    indices = dict(ndvi=ndvi, ndwi=ndwi, mndwi=mndwi, ndbi=ndbi, bsi=bsi)

    return RemovalResult(
        combined_mask   = combined,
        vegetation_mask = veg_mask,
        water_mask      = water_mask,
        buildup_mask    = build_mask,
        cloud_mask      = cloud_mask,
        snow_mask       = snow_mask,
        indices         = indices,
        stats           = stats,
    )


def apply_mask(image: np.ndarray,
               mask: np.ndarray,
               fill_value: float = np.nan) -> np.ndarray:
    """
    将掩膜应用于多波段影像。

    Parameters
    ----------
    image : np.ndarray  shape (bands, rows, cols)
    mask  : np.ndarray  shape (rows, cols)，True = 需剔除
    fill_value : float  填充值，默认 NaN

    Returns
    -------
    np.ndarray  与 image 同形状，受干扰像元已被 fill_value 替代
    """
    result = image.astype(np.float32).copy()
    result[:, mask] = fill_value
    return result


def print_stats(result: RemovalResult, sensor: str = "") -> None:
    """打印干扰剔除统计报告"""
    title = f"干扰剔除统计报告 [{sensor}]" if sensor else "干扰剔除统计报告"
    sep = "─" * 45
    print(f"\n{sep}")
    print(f"  {title}")
    print(sep)
    s = result.stats
    print(f"  总像元数         : {s['total_pixels']:,}")
    print(f"  植被覆盖         : {s['vegetation_pct']:6.2f} %")
    print(f"  水体覆盖         : {s['water_pct']:6.2f} %")
    print(f"  建筑/城镇        : {s['buildup_pct']:6.2f} %")
    print(f"  云覆盖           : {s['cloud_pct']:6.2f} %")
    print(f"  雪/冰覆盖        : {s['snow_pct']:6.2f} %")
    print(sep)
    print(f"  合计干扰像元     : {s['combined_pct']:6.2f} %")
    print(f"  有效地质像元     : {s['valid_geology_pct']:6.2f} %")
    print(sep)


# ─────────────────────────────────────────────
# 可选：保存掩膜 / 结果影像（依赖 rasterio）
# ─────────────────────────────────────────────

def save_masks(result: RemovalResult,
               output_dir: str,
               prefix: str = "mask") -> None:
    """
    将各类掩膜保存为 GeoTIFF（需安装 rasterio）。
    若 rasterio 未安装则保存为 .npy 文件。
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    masks = {
        "combined":   result.combined_mask,
        "vegetation": result.vegetation_mask,
        "water":      result.water_mask,
        "buildup":    result.buildup_mask,
        "cloud":      result.cloud_mask,
        "snow":       result.snow_mask,
    }

    try:
        import rasterio
        from rasterio.transform import from_bounds
        rows, cols = result.combined_mask.shape
        transform = from_bounds(0, 0, cols, rows, cols, rows)
        profile = dict(driver="GTiff", dtype="uint8",
                       count=1, height=rows, width=cols,
                       transform=transform, crs="EPSG:4326")
        for name, arr in masks.items():
            path = out / f"{prefix}_{name}.tif"
            with rasterio.open(path, "w", **profile) as dst:
                dst.write(arr.astype(np.uint8)[np.newaxis, :, :])
            print(f"  已保存: {path}")
    except ImportError:
        for name, arr in masks.items():
            path = out / f"{prefix}_{name}.npy"
            np.save(path, arr)
            print(f"  已保存: {path}")


# ─────────────────────────────────────────────
# 演示 / 单元测试
# ─────────────────────────────────────────────

def _demo():
    """用合成数据演示完整流程"""
    np.random.seed(42)
    bands, rows, cols = 7, 256, 256
    image = np.random.uniform(0.01, 0.35, (bands, rows, cols)).astype(np.float32)

    # 模拟典型地物反射率特征
    # 植被区块（高 NIR，低 Red）
    image[LANDSAT8_BANDS.nir,  50:100, 50:100]  = 0.45
    image[LANDSAT8_BANDS.red,  50:100, 50:100]  = 0.05

    # 水体区块（高 Green，低 NIR）
    image[LANDSAT8_BANDS.green, 120:160, 120:160] = 0.30
    image[LANDSAT8_BANDS.nir,   120:160, 120:160] = 0.02

    # 云区块（高 Blue + 高 NIR）
    image[LANDSAT8_BANDS.blue, 200:240, 200:240] = 0.40
    image[LANDSAT8_BANDS.nir,  200:240, 200:240] = 0.38

    result = remove_interference(image, band_cfg=LANDSAT8_BANDS)
    print_stats(result, sensor="Landsat 8 (合成数据演示)")

    # 应用掩膜，获得干净地质影像
    clean = apply_mask(image, result.combined_mask)
    valid_px = np.sum(~np.isnan(clean[0]))
    print(f"\n  剔除后有效像元数: {valid_px:,} / {rows * cols:,}")

    # 保存到当前目录
    save_masks(result, output_dir="./output_masks", prefix="geo")
    print("\n  演示完成。")


if __name__ == "__main__":
    _demo()
