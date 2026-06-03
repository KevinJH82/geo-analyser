"""
几何校正模块 - 对遥感影像进行地理配准和几何畸变校正
支持 Landsat 8/9、Sentinel-2 等多光谱遥感数据，可选依赖 rasterio 和 OpenCV
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Tuple, List
import warnings

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# 坐标系与投影配置
# ─────────────────────────────────────────────

@dataclass
class ProjectionConfig:
    """投影和坐标系配置"""
    crs: str = "EPSG:4326"          # 默认 WGS84
    pixel_size: float = 30.0        # 像元大小（米）
    upper_left_x: float = 0.0       # 左上角 X 坐标
    upper_left_y: float = 100.0     # 左上角 Y 坐标


LANDSAT8_PROJECTION = ProjectionConfig(
    crs="EPSG:4326",
    pixel_size=30.0,
    upper_left_x=0.0,
    upper_left_y=100.0
)

SENTINEL2_PROJECTION = ProjectionConfig(
    crs="EPSG:4326",
    pixel_size=10.0,
    upper_left_x=0.0,
    upper_left_y=100.0
)


# ─────────────────────────────────────────────
# 几何变换核心函数
# ─────────────────────────────────────────────

@dataclass
class GeoTransform:
    """地理变换参数（兼容 GDAL/rasterio 格式）"""
    origin_x: float              # 左上角 X
    pixel_width: float           # 像元宽度（X 方向）
    rotation_x: float = 0.0      # X 旋转（通常为 0）
    origin_y: float = 0.0        # 左上角 Y
    rotation_y: float = 0.0      # Y 旋转（通常为 0）
    pixel_height: float = -30.0  # 像元高度（Y 方向，通常为负）

    def to_tuple(self) -> Tuple[float, float, float, float, float, float]:
        """转换为 GDAL 标准 6 元组"""
        return (self.origin_x, self.pixel_width, self.rotation_x,
                self.origin_y, self.rotation_y, self.pixel_height)

    @classmethod
    def from_bounds(cls, ul_x: float, ul_y: float,
                   pixel_size: float = 30.0) -> "GeoTransform":
        """从边界和像元大小创建"""
        return cls(
            origin_x=ul_x,
            pixel_width=pixel_size,
            origin_y=ul_y,
            pixel_height=-pixel_size
        )


def pixel_to_geo(row: int, col: int,
                geo_transform: GeoTransform) -> Tuple[float, float]:
    """
    像素坐标转地理坐标（WGS84 或其他）

    Parameters
    ----------
    row, col : int
        像素行列号（从 0 开始）
    geo_transform : GeoTransform
        地理变换参数

    Returns
    -------
    (geo_x, geo_y) : Tuple[float, float]
        地理坐标
    """
    geo_x = geo_transform.origin_x + col * geo_transform.pixel_width
    geo_y = geo_transform.origin_y + row * geo_transform.pixel_height
    return geo_x, geo_y


def geo_to_pixel(geo_x: float, geo_y: float,
                geo_transform: GeoTransform) -> Tuple[int, int]:
    """
    地理坐标转像素坐标

    Parameters
    ----------
    geo_x, geo_y : float
        地理坐标
    geo_transform : GeoTransform
        地理变换参数

    Returns
    -------
    (row, col) : Tuple[int, int]
        像素坐标
    """
    col = int((geo_x - geo_transform.origin_x) / geo_transform.pixel_width)
    row = int((geo_y - geo_transform.origin_y) / geo_transform.pixel_height)
    return row, col


# ─────────────────────────────────────────────
# 基准点配准
# ─────────────────────────────────────────────

@dataclass
class GroundControlPoint:
    """地面控制点（GCP）"""
    pixel_x: float      # 影像像素 X
    pixel_y: float      # 影像像素 Y
    geo_x: float        # 已知地理坐标 X
    geo_y: float        # 已知地理坐标 Y


def compute_affine_transform(gcps: List[GroundControlPoint]) -> np.ndarray:
    """
    通过最小二乘法从 GCP 计算仿射变换矩阵（2x3）
    变换方程: [geo_x, geo_y]^T = A @ [pixel_x, pixel_y, 1]^T

    Parameters
    ----------
    gcps : List[GroundControlPoint]
        至少 3 个地面控制点

    Returns
    -------
    np.ndarray
        2×3 仿射变换矩阵
    """
    if len(gcps) < 3:
        raise ValueError("需要至少 3 个地面控制点")

    n = len(gcps)
    # 构造矩阵 [px, py, 1] 的堆叠
    A = np.column_stack([
        np.array([gcp.pixel_x for gcp in gcps]),
        np.array([gcp.pixel_y for gcp in gcps]),
        np.ones(n)
    ])

    # 地理坐标
    b_x = np.array([gcp.geo_x for gcp in gcps])
    b_y = np.array([gcp.geo_y for gcp in gcps])

    # 最小二乘求解
    coef_x, _, _, _ = np.linalg.lstsq(A, b_x, rcond=None)
    coef_y, _, _, _ = np.linalg.lstsq(A, b_y, rcond=None)

    # 返回 2×3 矩阵
    return np.vstack([coef_x, coef_y])


def apply_affine_transform(image: np.ndarray,
                          affine_matrix: np.ndarray,
                          output_shape: Optional[Tuple[int, int]] = None,
                          order: int = 1) -> np.ndarray:
    """
    使用仿射变换重采样影像（基于 scipy）

    Parameters
    ----------
    image : np.ndarray
        输入影像，shape (bands, rows, cols)
    affine_matrix : np.ndarray
        2×3 仿射变换矩阵
    output_shape : Tuple[int, int]
        输出影像大小 (rows, cols)，默认保持原大小
    order : int
        插值阶数（0=最近邻，1=双线性，3=三次）

    Returns
    -------
    np.ndarray
        变换后的影像
    """
    try:
        from scipy import ndimage
    except ImportError:
        raise ImportError("需要安装 scipy，执行: pip install scipy")

    if output_shape is None:
        output_shape = image.shape[1:]

    bands = image.shape[0]
    result = np.zeros((bands, *output_shape), dtype=image.dtype)

    # 构造逆变换矩阵（从目标像素反推源像素）
    try:
        inv_matrix = np.linalg.inv(
            np.vstack([affine_matrix, [0, 0, 1]])
        )[:2, :3]
    except np.linalg.LinAlgError:
        raise ValueError("仿射矩阵奇异，无法求逆")

    # 对每个波段应用变换
    for b in range(bands):
        coords = np.array([
            inv_matrix[0, 0] * np.arange(output_shape[1]) +
            inv_matrix[0, 1] * np.arange(output_shape[0])[:, np.newaxis] +
            inv_matrix[0, 2],
            inv_matrix[1, 0] * np.arange(output_shape[1]) +
            inv_matrix[1, 1] * np.arange(output_shape[0])[:, np.newaxis] +
            inv_matrix[1, 2]
        ])
        result[b] = ndimage.map_coordinates(
            image[b], coords, order=order, mode="constant", cval=0
        )

    return result.astype(image.dtype)


# ─────────────────────────────────────────────
# 多项式变换（高阶校正）
# ─────────────────────────────────────────────

def compute_polynomial_transform(gcps: List[GroundControlPoint],
                                degree: int = 2) -> dict:
    """
    通过最小二乘法计算多项式变换（2-3 阶）
    用于非线性几何畸变（如镜头畸变、轨道倾斜）

    Parameters
    ----------
    gcps : List[GroundControlPoint]
        地面控制点列表
    degree : int
        多项式阶数（1=线性/仿射，2=二次，3=三次）

    Returns
    -------
    dict
        包含 X 和 Y 方向多项式系数的字典
    """
    if len(gcps) < (degree + 1) ** 2:
        raise ValueError(
            f"度数 {degree} 需要至少 {(degree + 1) ** 2} 个 GCP"
        )

    px = np.array([gcp.pixel_x for gcp in gcps])
    py = np.array([gcp.pixel_y for gcp in gcps])
    gx = np.array([gcp.geo_x for gcp in gcps])
    gy = np.array([gcp.geo_y for gcp in gcps])

    # 构造多项式特征矩阵
    features = []
    for d in range(degree + 1):
        for i in range(d + 1):
            features.append(px ** (d - i) * py ** i)

    A = np.column_stack(features)

    # 最小二乘求解
    coef_x, _, _, _ = np.linalg.lstsq(A, gx, rcond=None)
    coef_y, _, _, _ = np.linalg.lstsq(A, gy, rcond=None)

    return {
        "degree": degree,
        "coef_x": coef_x,
        "coef_y": coef_y,
    }


def apply_polynomial_transform(image: np.ndarray,
                              poly_params: dict,
                              output_shape: Optional[Tuple[int, int]] = None,
                              order: int = 1) -> np.ndarray:
    """
    应用多项式变换重采样影像

    Parameters
    ----------
    image : np.ndarray
        输入影像
    poly_params : dict
        多项式参数（来自 compute_polynomial_transform）
    output_shape : Tuple[int, int]
        输出大小
    order : int
        插值阶数

    Returns
    -------
    np.ndarray
        变换后的影像
    """
    try:
        from scipy import ndimage
    except ImportError:
        raise ImportError("需要安装 scipy")

    if output_shape is None:
        output_shape = image.shape[1:]

    degree = poly_params["degree"]
    coef_x = poly_params["coef_x"]
    coef_y = poly_params["coef_y"]

    # 生成目标网格
    row_grid, col_grid = np.mgrid[0:output_shape[0], 0:output_shape[1]]

    # 计算多项式变换
    src_x = np.zeros_like(row_grid, dtype=np.float32)
    src_y = np.zeros_like(row_grid, dtype=np.float32)

    coef_idx = 0
    for d in range(degree + 1):
        for i in range(d + 1):
            term = col_grid ** (d - i) * row_grid ** i
            src_x += coef_x[coef_idx] * term
            src_y += coef_y[coef_idx] * term
            coef_idx += 1

    bands = image.shape[0]
    result = np.zeros((bands, *output_shape), dtype=image.dtype)

    for b in range(bands):
        result[b] = ndimage.map_coordinates(
            image[b],
            [src_y, src_x],
            order=order,
            mode="constant",
            cval=0
        )

    return result.astype(image.dtype)


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

@dataclass
class CorrectionResult:
    """几何校正结果"""
    corrected_image: np.ndarray        # 校正后的影像
    geo_transform: GeoTransform        # 地理变换参数
    affine_matrix: Optional[np.ndarray] = None  # 仿射矩阵（若使用）
    poly_params: Optional[dict] = None          # 多项式参数（若使用）
    gcps_used: Optional[List[GroundControlPoint]] = None
    stats: dict = field(default_factory=dict)


def geometric_correction(
    image: np.ndarray,
    gcps: List[GroundControlPoint],
    projection: ProjectionConfig = LANDSAT8_PROJECTION,
    method: str = "affine",
    poly_degree: int = 2,
    output_shape: Optional[Tuple[int, int]] = None,
    interpolation: str = "bilinear",
) -> CorrectionResult:
    """
    对多光谱影像执行几何校正。

    Parameters
    ----------
    image : np.ndarray
        输入影像，shape (bands, rows, cols)，值需归一化到 [0, 1]
    gcps : List[GroundControlPoint]
        地面控制点列表（至少 3 个）
    projection : ProjectionConfig
        投影配置
    method : str
        校正方法："affine" 或 "polynomial"
    poly_degree : int
        多项式度数（仅当 method="polynomial" 时有效）
    output_shape : Tuple[int, int]
        输出影像大小，默认保持原大小
    interpolation : str
        插值方法："nearest"、"bilinear" 或 "cubic"

    Returns
    -------
    CorrectionResult
        校正后的影像和变换参数
    """
    if output_shape is None:
        output_shape = image.shape[1:]

    # 插值阶数映射
    order_map = {"nearest": 0, "bilinear": 1, "cubic": 3}
    order = order_map.get(interpolation, 1)

    if method == "affine":
        affine_matrix = compute_affine_transform(gcps)
        corrected = apply_affine_transform(
            image, affine_matrix, output_shape, order
        )
        poly_params = None
    elif method == "polynomial":
        poly_params = compute_polynomial_transform(gcps, poly_degree)
        corrected = apply_polynomial_transform(
            image, poly_params, output_shape, order
        )
        affine_matrix = None
    else:
        raise ValueError("method 必须是 'affine' 或 'polynomial'")

    # 构造地理变换
    geo_transform = GeoTransform.from_bounds(
        projection.upper_left_x,
        projection.upper_left_y,
        projection.pixel_size
    )

    # 统计信息
    stats = {
        "method": method,
        "gcps_count": len(gcps),
        "output_shape": output_shape,
        "projection": projection.crs,
        "pixel_size": projection.pixel_size,
    }

    return CorrectionResult(
        corrected_image=corrected,
        geo_transform=geo_transform,
        affine_matrix=affine_matrix,
        poly_params=poly_params,
        gcps_used=gcps,
        stats=stats,
    )


def print_correction_stats(result: CorrectionResult) -> None:
    """打印几何校正统计报告"""
    sep = "─" * 50
    print(f"\n{sep}")
    print(f"  几何校正统计报告")
    print(sep)
    s = result.stats
    print(f"  校正方法         : {s['method']}")
    print(f"  使用 GCP 数       : {s['gcps_count']}")
    print(f"  输出影像大小     : {s['output_shape']}")
    print(f"  投影系统         : {s['projection']}")
    print(f"  像元大小         : {s['pixel_size']} m")
    print(sep)

    gt = result.geo_transform
    print(f"  地理变换参数:")
    print(f"    左上角 X       : {gt.origin_x:.6f}")
    print(f"    左上角 Y       : {gt.origin_y:.6f}")
    print(f"    像元宽度       : {gt.pixel_width:.6f}")
    print(f"    像元高度       : {gt.pixel_height:.6f}")
    print(sep)


# ─────────────────────────────────────────────
# 可选：保存校正结果（依赖 rasterio）
# ─────────────────────────────────────────────

def save_corrected_image(result: CorrectionResult,
                        output_path: str,
                        save_geotiff: bool = True) -> None:
    """
    保存几何校正后的影像为 GeoTIFF（需 rasterio）
    或 .npy 格式。

    Parameters
    ----------
    result : CorrectionResult
        校正结果
    output_path : str
        输出路径
    save_geotiff : bool
        是否尝试保存为 GeoTIFF
    """
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if save_geotiff:
        try:
            import rasterio
            from rasterio.transform import Affine

            bands, rows, cols = result.corrected_image.shape
            gt = result.geo_transform
            transform = Affine(
                gt.pixel_width, gt.rotation_x, gt.origin_x,
                gt.rotation_y, gt.pixel_height, gt.origin_y
            )

            profile = dict(
                driver="GTiff",
                dtype=result.corrected_image.dtype,
                count=bands,
                height=rows,
                width=cols,
                transform=transform,
                crs="EPSG:4326"
            )

            with rasterio.open(out_path, "w", **profile) as dst:
                for b in range(bands):
                    dst.write(result.corrected_image[b], b + 1)

            print(f"  已保存 GeoTIFF: {out_path}")
        except ImportError:
            print("  rasterio 未安装，保存为 .npy 格式")
            np.save(out_path.with_suffix(".npy"), result.corrected_image)
    else:
        np.save(out_path.with_suffix(".npy"), result.corrected_image)
        print(f"  已保存 .npy: {out_path.with_suffix('.npy')}")


# ─────────────────────────────────────────────
# 演示 / 单元测试
# ─────────────────────────────────────────────

def _demo():
    """用合成数据演示完整的几何校正流程"""
    np.random.seed(42)
    bands, rows, cols = 7, 256, 256

    # 生成合成影像
    image = np.random.uniform(0.01, 0.35, (bands, rows, cols)).astype(np.float32)

    # 模拟地物特征
    image[:, 50:100, 50:100] = 0.25   # 特征区域 1
    image[:, 150:200, 150:200] = 0.15 # 特征区域 2

    # 创建 GCP（已知的像素-地理坐标对应关系）
    # 对于仿射变换，3 个 GCP 即可
    gcps_affine = [
        GroundControlPoint(pixel_x=0, pixel_y=0, geo_x=0.0, geo_y=100.0),
        GroundControlPoint(pixel_x=cols, pixel_y=0, geo_x=cols*30/111000, geo_y=100.0),
        GroundControlPoint(pixel_x=0, pixel_y=rows, geo_x=0.0, geo_y=100-rows*30/111000),
    ]

    # 对于二阶多项式变换，需要至少 9 个 GCP
    gcps_poly = [
        GroundControlPoint(pixel_x=0, pixel_y=0, geo_x=0.0, geo_y=100.0),
        GroundControlPoint(pixel_x=cols//2, pixel_y=0, geo_x=(cols//2)*30/111000, geo_y=100.0),
        GroundControlPoint(pixel_x=cols, pixel_y=0, geo_x=cols*30/111000, geo_y=100.0),
        GroundControlPoint(pixel_x=0, pixel_y=rows//2, geo_x=0.0, geo_y=100-(rows//2)*30/111000),
        GroundControlPoint(pixel_x=cols//2, pixel_y=rows//2,
                          geo_x=(cols//2)*30/111000, geo_y=100-(rows//2)*30/111000),
        GroundControlPoint(pixel_x=cols, pixel_y=rows//2,
                          geo_x=cols*30/111000, geo_y=100-(rows//2)*30/111000),
        GroundControlPoint(pixel_x=0, pixel_y=rows, geo_x=0.0, geo_y=100-rows*30/111000),
        GroundControlPoint(pixel_x=cols//2, pixel_y=rows,
                          geo_x=(cols//2)*30/111000, geo_y=100-rows*30/111000),
        GroundControlPoint(pixel_x=cols, pixel_y=rows,
                          geo_x=cols*30/111000, geo_y=100-rows*30/111000),
    ]

    # 执行仿射校正
    print("\n" + "="*50)
    print("  仿射变换校正")
    print("="*50)
    result_affine = geometric_correction(
        image, gcps_affine, method="affine", interpolation="bilinear"
    )
    print_correction_stats(result_affine)

    # 执行多项式校正
    print("\n" + "="*50)
    print("  二阶多项式校正")
    print("="*50)
    result_poly = geometric_correction(
        image, gcps_poly, method="polynomial", poly_degree=2, interpolation="bilinear"
    )
    print_correction_stats(result_poly)

    # 保存结果
    save_corrected_image(result_affine, "./output_masks/corrected_affine.tif")
    save_corrected_image(result_poly, "./output_masks/corrected_poly.tif")

    print("\n  演示完成。")


if __name__ == "__main__":
    _demo()
