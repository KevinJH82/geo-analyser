"""
大气校正模块 - 对遥感影像进行大气校正，去除大气干扰影响
支持 Landsat 8/9、Sentinel-2 等多光谱遥感数据
主要方法：DOS（暗目标消减）、相对辐射校正、大气透过率计算
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict
import warnings

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# 辐射标定参数（Landsat 8/9）
# ─────────────────────────────────────────────

@dataclass
class RadianceCoefficients:
    """DN 值转辐射亮度的线性标定系数（Landsat 8/9）"""
    ml: float  # 辐射倍数因子（Radiance Mult Band）
    al: float  # 辐射加法因子（Radiance Add Band）
    band_name: str = "unknown"


# Landsat 8/9 辐射标定参数（示例值，实际应从元数据读取）
LANDSAT8_RADIANCE = {
    1: RadianceCoefficients(ml=0.0003342, al=0.1, band_name="Coastal"),
    2: RadianceCoefficients(ml=0.0003344, al=0.1, band_name="Blue"),
    3: RadianceCoefficients(ml=0.0003303, al=0.1, band_name="Green"),
    4: RadianceCoefficients(ml=0.0002699, al=0.1, band_name="Red"),
    5: RadianceCoefficients(ml=0.0002837, al=0.1, band_name="NIR"),
    6: RadianceCoefficients(ml=0.0000515, al=0.1, band_name="SWIR1"),
    7: RadianceCoefficients(ml=0.0000705, al=0.1, band_name="SWIR2"),
}

LANDSAT8_REFLECTANCE = {
    1: RadianceCoefficients(ml=0.0002, al=0.0, band_name="Coastal"),
    2: RadianceCoefficients(ml=0.0002, al=0.0, band_name="Blue"),
    3: RadianceCoefficients(ml=0.0002, al=0.0, band_name="Green"),
    4: RadianceCoefficients(ml=0.0002, al=0.0, band_name="Red"),
    5: RadianceCoefficients(ml=0.0001, al=0.0, band_name="NIR"),
    6: RadianceCoefficients(ml=0.00005, al=0.0, band_name="SWIR1"),
    7: RadianceCoefficients(ml=0.00005, al=0.0, band_name="SWIR2"),
}


# ─────────────────────────────────────────────
# 波段配置与太阳几何参数
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


@dataclass
class SolarGeometry:
    """太阳几何参数"""
    zenith: float      # 太阳天顶角（度）
    azimuth: float     # 太阳方位角（度）
    elevation: float   # 太阳高度角（度）= 90 - zenith
    distance: float    # 日地距修正因子（相对标准距离）

    @classmethod
    def from_zenith(cls, zenith: float) -> "SolarGeometry":
        """从太阳天顶角创建"""
        return cls(
            zenith=zenith,
            azimuth=0.0,
            elevation=90.0 - zenith,
            distance=1.0
        )


LANDSAT8_BANDS = BandConfig(
    blue=1, green=2, red=3, nir=4, swir1=5, swir2=6, name="Landsat8/9"
)

SENTINEL2_BANDS = BandConfig(
    blue=1, green=2, red=3, nir=7, swir1=9, swir2=10, name="Sentinel-2"
    # 读入顺序: B01(0),B02(1),B03(2),B04(3),B05(4),B06(5),B07(6),B08(7),B8A(8),B11(9),B12(10)
)


# ─────────────────────────────────────────────
# 辐射标定函数
# ─────────────────────────────────────────────

def dn_to_radiance(dn: np.ndarray,
                  ml: float, al: float) -> np.ndarray:
    """
    DN 值转辐射亮度（Landsat 8/9）
    L_λ = ml * DN + al

    Parameters
    ----------
    dn : np.ndarray
        DN 值（整型，范围 [0, 65535]）
    ml : float
        辐射倍数因子
    al : float
        辐射加法因子

    Returns
    -------
    np.ndarray
        辐射亮度值
    """
    return ml * dn.astype(np.float32) + al


def radiance_to_reflectance(radiance: np.ndarray,
                           solar_zenith: float,
                           distance: float = 1.0,
                           esun: float = 1500.0) -> np.ndarray:
    """
    辐射亮度转反射率（TOA 反射率）
    ρ_TOA = π * L_λ / (ESUN_λ * cos(θ_s) * d^2)

    Parameters
    ----------
    radiance : np.ndarray
        辐射亮度值
    solar_zenith : float
        太阳天顶角（度）
    distance : float
        日地距修正因子
    esun : float
        太阳常数各波段值

    Returns
    -------
    np.ndarray
        TOA 反射率，范围 [0, 1]
    """
    cos_zenith = np.cos(np.radians(solar_zenith))
    rho_toa = (np.pi * radiance) / (esun * cos_zenith * distance ** 2)
    return np.clip(rho_toa.astype(np.float32), 0.0, 1.0)


# ─────────────────────────────────────────────
# DOS 大气校正（暗目标消减）
# ─────────────────────────────────────────────

@dataclass
class DOSParameters:
    """DOS 校正参数"""
    ndvi_threshold: float = 0.4   # NDVI 阈值识别暗目标
    percentile: float = 0.01       # 暗目标 DN 百分位（通常 0.01-0.1%）
    sky_radiance: Optional[float] = None  # 大气顶部辐射值


def estimate_dos_correction(
    image: np.ndarray,
    band_cfg: BandConfig,
    band_idx: int,
    solar_zenith: float,
    esun_dict: Dict[int, float],
    params: Optional[DOSParameters] = None,
) -> Tuple[np.ndarray, Dict]:
    """
    基于 DOS 方法估计大气校正

    Parameters
    ----------
    image : np.ndarray
        输入影像，shape (bands, rows, cols)，DN 值
    band_cfg : BandConfig
        波段配置
    band_idx : int
        要校正的波段索引
    solar_zenith : float
        太阳天顶角（度）
    esun_dict : Dict[int, float]
        各波段太阳常数
    params : DOSParameters
        DOS 参数

    Returns
    -------
    (corrected, stats) : Tuple[np.ndarray, Dict]
        校正后的反射率和统计信息
    """
    if params is None:
        params = DOSParameters()

    band = image[band_idx].astype(np.float32)

    # 计算 NDVI 识别植被/暗目标区域
    def safe_divide(a, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(b != 0, a / b, 0.0)

    nir = image[band_cfg.nir].astype(np.float32)
    red = image[band_cfg.red].astype(np.float32)
    ndvi = safe_divide(nir - red, nir + red)

    # 识别暗目标（低 NDVI 且 DN 值低）
    dark_mask = ndvi < -0.1  # 非植被区域
    dark_pixels = band[dark_mask]

    if len(dark_pixels) == 0:
        # 如果没有暗目标，取全局最小百分位
        min_dn = np.percentile(band, params.percentile * 100)
    else:
        min_dn = np.percentile(dark_pixels, params.percentile * 100)

    # 转换为辐射亮度
    ml = LANDSAT8_RADIANCE[band_idx + 1].ml
    al = LANDSAT8_RADIANCE[band_idx + 1].al
    min_radiance = dn_to_radiance(np.array([min_dn]), ml, al)[0]

    # 转换为反射率（TOA 反射率）
    esun = esun_dict.get(band_idx + 1, 1500.0)
    min_reflectance = radiance_to_reflectance(
        np.array([min_radiance]), solar_zenith, esun=esun
    )[0]

    # DOS 校正
    corrected = band.copy()
    corrected = radiance_to_reflectance(
        dn_to_radiance(corrected, ml, al),
        solar_zenith,
        esun=esun
    )
    # 减去大气顶部辐射
    corrected = np.maximum(corrected - min_reflectance, 0.0)

    stats = {
        "band_idx": band_idx,
        "min_dn": float(min_dn),
        "min_reflectance": float(min_reflectance),
        "dark_pixels_count": int(len(dark_pixels)),
    }

    return corrected.astype(np.float32), stats


# ─────────────────────────────────────────────
# 相对辐射校正（相对 DN 值校正）
# ─────────────────────────────────────────────

def relative_radiometric_correction(
    image: np.ndarray,
    reference_band: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    相对辐射校正 - 消除条纹噪声和传感器间差异
    基于每条扫描线的平均 DN 值进行归一化

    Parameters
    ----------
    image : np.ndarray
        输入影像，shape (bands, rows, cols)
    reference_band : np.ndarray
        参考波段用于计算扫描线均值

    Returns
    -------
    np.ndarray
        校正后的影像
    """
    if reference_band is None:
        reference_band = image[1]  # 使用绿波段作为参考

    bands, rows, cols = image.shape
    corrected = image.astype(np.float32).copy()

    # 计算每行的平均 DN
    line_means = np.nanmean(reference_band.astype(np.float32), axis=1)
    global_mean = np.nanmean(line_means)

    # 对每行进行相对校正
    for b in range(bands):
        for r in range(rows):
            if line_means[r] != 0:
                corrected[b, r, :] *= (global_mean / line_means[r])

    return corrected.astype(image.dtype)


# ─────────────────────────────────────────────
# 简化大气校正（相对反射率）
# ─────────────────────────────────────────────

def simplified_atmospheric_correction(
    image: np.ndarray,
    band_cfg: BandConfig,
    solar_zenith: float = 30.0,
    distance: float = 1.0,
    esun_dict: Optional[Dict[int, float]] = None,
) -> np.ndarray:
    """
    简化大气校正 - 直接 DN→反射率转换
    适用于需要快速处理、精度要求不高的场景

    Parameters
    ----------
    image : np.ndarray
        输入影像（DN 值）
    band_cfg : BandConfig
        波段配置
    solar_zenith : float
        太阳天顶角（度）
    distance : float
        日地距修正因子
    esun_dict : Dict[int, float]
        各波段太阳常数

    Returns
    -------
    np.ndarray
        TOA 反射率，范围 [0, 1]
    """
    if esun_dict is None:
        # Landsat 8/9 TOA 反射率标定系数（简化）
        esun_dict = {
            1: 1895, 2: 1941, 3: 1822, 4: 1533,
            5: 1039, 6: 374.8, 7: 224.4
        }

    bands = image.shape[0]
    corrected = np.zeros_like(image, dtype=np.float32)

    cos_zenith = np.cos(np.radians(solar_zenith))

    for b in range(bands):
        band_dn = image[b].astype(np.float32)
        ml = LANDSAT8_RADIANCE[b + 1].ml
        al = LANDSAT8_RADIANCE[b + 1].al

        # DN → 辐射亮度
        radiance = dn_to_radiance(band_dn, ml, al)

        # 辐射亮度 → 反射率
        esun = esun_dict.get(b + 1, 1500.0)
        reflectance = (np.pi * radiance) / (esun * cos_zenith * distance ** 2)
        corrected[b] = np.clip(reflectance, 0.0, 1.0)

    return corrected


# ─────────────────────────────────────────────
# 大气透过率与气溶胶光学厚度
# ─────────────────────────────────────────────

def calculate_atmospheric_transmission(
    aot550: float = 0.1,
    altitude: float = 0.0,
    zenith_angle: float = 30.0,
) -> float:
    """
    计算大气透过率（简化模型）
    基于气溶胶光学厚度（AOT）

    Parameters
    ----------
    aot550 : float
        550nm 处的气溶胶光学厚度
    altitude : float
        地表高度（km）
    zenith_angle : float
        天顶角（度）

    Returns
    -------
    float
        大气透过率 [0, 1]
    """
    # 简化的 Kasten-Czeplak 模型
    air_mass = 1.0 / (np.cos(np.radians(zenith_angle)) + 0.50572 * (96.07995 - zenith_angle) ** -1.6364)

    # 考虑高度的气溶胶光学厚度
    aot_adj = aot550 * np.exp(-altitude / 8.0)

    # 透过率
    transmission = np.exp(-aot_adj * air_mass)
    return float(np.clip(transmission, 0.0, 1.0))


def estimate_surface_reflectance(
    toa_reflectance: np.ndarray,
    transmission: float,
    sky_radiance: float = 0.05,
) -> np.ndarray:
    """
    从 TOA 反射率估计地表反射率

    Parameters
    ----------
    toa_reflectance : np.ndarray
        TOA 反射率
    transmission : float
        大气透过率
    sky_radiance : float
        漫反射辐射比例

    Returns
    -------
    np.ndarray
        地表反射率
    """
    if transmission <= 0:
        transmission = 0.01

    surface = (toa_reflectance - sky_radiance) / transmission
    return np.clip(surface.astype(np.float32), 0.0, 1.0)


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

@dataclass
class CorrectionResult:
    """大气校正结果"""
    surface_reflectance: np.ndarray    # 地表反射率
    toa_reflectance: np.ndarray        # TOA 反射率
    method: str                        # 校正方法名称
    solar_geometry: SolarGeometry      # 太阳几何参数
    stats: dict = field(default_factory=dict)  # 统计信息


def atmospheric_correction(
    image: np.ndarray,
    band_cfg: BandConfig = LANDSAT8_BANDS,
    solar_zenith: float = 30.0,
    distance: float = 1.0,
    method: str = "dos",
    dos_params: Optional[DOSParameters] = None,
    aot550: float = 0.1,
    altitude: float = 0.0,
    esun_dict: Optional[Dict[int, float]] = None,
) -> CorrectionResult:
    """
    对多光谱影像执行完整的大气校正流程。

    Parameters
    ----------
    image : np.ndarray
        输入影像，shape (bands, rows, cols)，DN 值
    band_cfg : BandConfig
        波段索引配置
    solar_zenith : float
        太阳天顶角（度）
    distance : float
        日地距修正因子
    method : str
        校正方法："dos"（暗目标消减）或 "simple"（简化）
    dos_params : DOSParameters
        DOS 参数
    aot550 : float
        550nm 处的气溶胶光学厚度
    altitude : float
        地表高度（km）
    esun_dict : Dict[int, float]
        各波段太阳常数

    Returns
    -------
    CorrectionResult
        校正结果
    """
    if esun_dict is None:
        esun_dict = {
            1: 1895, 2: 1941, 3: 1822, 4: 1533,
            5: 1039, 6: 374.8, 7: 224.4
        }

    solar_geom = SolarGeometry.from_zenith(solar_zenith)

    # 先计算 TOA 反射率
    toa_reflectance = simplified_atmospheric_correction(
        image, band_cfg, solar_zenith, distance, esun_dict
    )

    # 根据选择的方法进行大气校正
    if method == "dos":
        # DOS 方法校正各波段
        surface_reflectance = np.zeros_like(toa_reflectance)
        stats_list = []

        for b in range(image.shape[0]):
            corrected, stats = estimate_dos_correction(
                image, band_cfg, b, solar_zenith, esun_dict, dos_params
            )
            surface_reflectance[b] = corrected
            stats_list.append(stats)

        stats = {
            "method": "dos",
            "solar_zenith": solar_zenith,
            "distance": distance,
            "band_stats": stats_list,
        }

    elif method == "simple":
        # 简化方法：仅使用透过率
        transmission = calculate_atmospheric_transmission(
            aot550, altitude, solar_zenith
        )
        surface_reflectance = estimate_surface_reflectance(
            toa_reflectance, transmission
        )

        stats = {
            "method": "simple",
            "solar_zenith": solar_zenith,
            "distance": distance,
            "aot550": aot550,
            "transmission": transmission,
        }

    else:
        raise ValueError("method 必须是 'dos' 或 'simple'")

    return CorrectionResult(
        surface_reflectance=surface_reflectance,
        toa_reflectance=toa_reflectance,
        method=method,
        solar_geometry=solar_geom,
        stats=stats,
    )


def print_correction_stats(result: CorrectionResult, sensor: str = "") -> None:
    """打印大气校正统计报告"""
    title = f"大气校正统计报告 [{sensor}]" if sensor else "大气校正统计报告"
    sep = "─" * 50
    print(f"\n{sep}")
    print(f"  {title}")
    print(sep)

    s = result.stats
    print(f"  校正方法         : {s['method'].upper()}")
    print(f"  太阳天顶角       : {s['solar_zenith']:.2f}°")
    print(f"  日地距修正因子   : {s['distance']:.4f}")

    if s["method"] == "dos":
        print(f"  处理波段数       : {len(s['band_stats'])}")
        print(f"  暗目标检测:")
        for stat in s["band_stats"]:
            print(f"    波段 {stat['band_idx'] + 1}: "
                  f"min_DN={stat['min_dn']:.1f}, "
                  f"dark_pixels={stat['dark_pixels_count']:,}")
    elif s["method"] == "simple":
        print(f"  气溶胶光学厚度   : {s['aot550']:.4f}")
        print(f"  大气透过率       : {s['transmission']:.4f}")

    print(sep)


# ─────────────────────────────────────────────
# 可选：保存校正结果（依赖 rasterio）
# ─────────────────────────────────────────────

def save_corrected_image(result: CorrectionResult,
                        output_path: str,
                        output_type: str = "surface",
                        save_geotiff: bool = True) -> None:
    """
    保存大气校正后的影像为 GeoTIFF 或 .npy 格式。

    Parameters
    ----------
    result : CorrectionResult
        校正结果
    output_path : str
        输出路径
    output_type : str
        输出类型："surface"（地表反射率）或 "toa"（TOA 反射率）
    save_geotiff : bool
        是否尝试保存为 GeoTIFF
    """
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    image = result.surface_reflectance if output_type == "surface" else result.toa_reflectance
    image_uint8 = (image * 255).astype(np.uint8)

    if save_geotiff:
        try:
            import rasterio
            from rasterio.transform import from_bounds

            bands, rows, cols = image.shape
            transform = from_bounds(0, 0, cols, rows, cols, rows)
            profile = dict(
                driver="GTiff",
                dtype="uint8",
                count=bands,
                height=rows,
                width=cols,
                transform=transform,
                crs="EPSG:4326"
            )

            with rasterio.open(out_path, "w", **profile) as dst:
                for b in range(bands):
                    dst.write(image_uint8[b], b + 1)

            print(f"  已保存 GeoTIFF: {out_path}")
        except ImportError:
            print("  rasterio 未安装，保存为 .npy 格式")
            np.save(out_path.with_suffix(".npy"), image)
    else:
        np.save(out_path.with_suffix(".npy"), image)
        print(f"  已保存 .npy: {out_path.with_suffix('.npy')}")


# ─────────────────────────────────────────────
# 演示 / 单元测试
# ─────────────────────────────────────────────

def _demo():
    """用合成数据演示完整的大气校正流程"""
    np.random.seed(42)
    bands, rows, cols = 7, 256, 256

    # 生成合成 DN 值影像（模拟 Landsat 8/9）
    # DN 值通常范围 [0, 65535]
    image = np.random.uniform(100, 5000, (bands, rows, cols)).astype(np.uint16)

    # 模拟典型地物反射率特征
    # 植被区块
    image[LANDSAT8_BANDS.nir,  50:100, 50:100]  = 4500
    image[LANDSAT8_BANDS.red,  50:100, 50:100]  = 1000

    # 水体区块
    image[LANDSAT8_BANDS.green, 120:160, 120:160] = 2000
    image[LANDSAT8_BANDS.nir,   120:160, 120:160] = 500

    # 建筑区块
    image[LANDSAT8_BANDS.red, 200:240, 200:240] = 3000
    image[LANDSAT8_BANDS.nir, 200:240, 200:240] = 2000

    # 执行 DOS 大气校正
    print("\n" + "="*50)
    print("  DOS 大气校正")
    print("="*50)
    result_dos = atmospheric_correction(
        image,
        band_cfg=LANDSAT8_BANDS,
        solar_zenith=30.0,
        distance=1.0,
        method="dos"
    )
    print_correction_stats(result_dos, sensor="Landsat 8 (合成数据)")

    # 执行简化大气校正
    print("\n" + "="*50)
    print("  简化大气校正")
    print("="*50)
    result_simple = atmospheric_correction(
        image,
        band_cfg=LANDSAT8_BANDS,
        solar_zenith=30.0,
        distance=1.0,
        method="simple",
        aot550=0.1
    )
    print_correction_stats(result_simple, sensor="Landsat 8 (合成数据)")

    # 输出反射率统计
    print(f"\n  DOS 方法结果:")
    print(f"    TOA 反射率范围   : [{result_dos.toa_reflectance.min():.4f}, "
          f"{result_dos.toa_reflectance.max():.4f}]")
    print(f"    地表反射率范围   : [{result_dos.surface_reflectance.min():.4f}, "
          f"{result_dos.surface_reflectance.max():.4f}]")

    print(f"\n  简化方法结果:")
    print(f"    TOA 反射率范围   : [{result_simple.toa_reflectance.min():.4f}, "
          f"{result_simple.toa_reflectance.max():.4f}]")
    print(f"    地表反射率范围   : [{result_simple.surface_reflectance.min():.4f}, "
          f"{result_simple.surface_reflectance.max():.4f}]")

    # 保存结果
    save_corrected_image(result_dos, "./output_masks/corrected_dos.tif")
    save_corrected_image(result_simple, "./output_masks/corrected_simple.tif")

    print("\n  演示完成。")


if __name__ == "__main__":
    _demo()
