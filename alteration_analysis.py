"""
蚀变异常检测算法层

两种方法:
  - 波段比值法 (ratio): 解析 alteration_deposit_db.json 中形如 "B5/B6"、"(B5+B7)/B6" 的表达式
  - Crosta PCA 法 (pca): 按 input_bands + pc_criterion(正/负载荷)自动识别目标 PC

支持传感器: Landsat 8/9、Sentinel-2、ASTER SWIR
所有阈值统计限定在 ROI 内进行。
"""

from __future__ import annotations

import ast
import re
import numpy as np
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# 兼容旧 API
from alteration_db import (
    get_recommended_targets,
    get_deposit_type_meta,
    normalize_sensor,
)
# 丛林模式:红边/地植物学胁迫指数 — 从 commons 共享层引入(纯函数,无循环依赖)
import os as _os, sys as _sys
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _REPO_ROOT not in _sys.path:
    _sys.path.insert(0, _REPO_ROOT)
from commons.spectral_indices import calc_ndvi, calc_ndre, calc_cire, calc_rep

# ─────────────────────────────────────────────
# 传感器波段配置 — JSON 中的 BN 编号 → image 数组通道索引
# ─────────────────────────────────────────────

# Landsat 8/9: JSON 用 B2(Blue) B5(NIR) B6(SWIR1) B7(SWIR2) 等
# 实际 image 是 7 波段顺序: B1 Coastal, B2 Blue, B3 Green, B4 Red, B5 NIR, B6 SWIR1, B7 SWIR2
LANDSAT8_BN_TO_IDX = {f"B{i}": i - 1 for i in range(1, 8)}

# Sentinel-2 L2A: 实际加载顺序 (B2,B3,B4,B5,B6,B7,B8,B8A,B9,B11,B12)
# JSON 中用 B2..B12 / B8A
SENTINEL2_BN_TO_IDX = {
    "B2": 0, "B3": 1, "B4": 2, "B5": 3, "B6": 4, "B7": 5,
    "B8": 6, "B8A": 7, "B9": 8, "B11": 9, "B12": 10,
}

# ASTER SWIR (6 波段加载): B4..B9 → 索引 0..5
ASTER_BN_TO_IDX = {f"B{i}": i - 4 for i in range(4, 10)}

# 兼容旧代码命名
LANDSAT_BANDS = {
    "coastal": 0, "blue": 1, "green": 2, "red": 3,
    "nir": 4, "swir1": 5, "swir2": 6,
}
SENTINEL2_BANDS = {
    "blue": 0, "green": 1, "red": 2,
    "nir": 6, "swir1": 9, "swir2": 10,
}
ASTER_SWIR_BANDS = {
    "b4": 0, "b5": 1, "b6": 2, "b7": 3, "b8": 4, "b9": 5,
}


def _bn_map(sensor_key: str) -> Dict[str, int]:
    if sensor_key == "Landsat8":
        return LANDSAT8_BN_TO_IDX
    if sensor_key == "Sentinel2":
        return SENTINEL2_BN_TO_IDX
    if sensor_key == "ASTER":
        return ASTER_BN_TO_IDX
    raise ValueError(f"不支持的传感器: {sensor_key}")


# ─────────────────────────────────────────────
# 安全表达式求值 (波段比值法)
# ─────────────────────────────────────────────

_BAND_TOKEN = re.compile(r"B\d+A?")  # 兼容 Sentinel-2 的 B8A

# AST 节点白名单
_ALLOWED_NODES = (
    ast.Expression, ast.BinOp, ast.UnaryOp, ast.Name, ast.Load,
    ast.Constant,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.USub, ast.UAdd,
)


def _validate_ast(node: ast.AST) -> None:
    for child in ast.walk(node):
        if not isinstance(child, _ALLOWED_NODES):
            raise ValueError(f"表达式包含不允许的节点: {type(child).__name__}")
        if isinstance(child, ast.Name) and not _BAND_TOKEN.fullmatch(child.id):
            raise ValueError(f"非法波段标识符: {child.id}")


def _safe_ratio(a: np.ndarray, b: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    denom = b.astype(np.float32)
    return np.where(np.abs(denom) > eps, a.astype(np.float32) / denom, np.nan).astype(np.float32)


class _SafeBinOpEval(ast.NodeVisitor):
    """白名单 AST 求值器。除法走 _safe_ratio,其他用 numpy 算子。"""

    def __init__(self, image: np.ndarray, sensor_key: str, bn_map: Optional[Dict[str, int]] = None):
        self.image = image
        self.bn_map = bn_map if bn_map is not None else _bn_map(sensor_key)
        self.sensor_key = sensor_key

    def visit_Expression(self, node):
        return self.visit(node.body)

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        if isinstance(node.op, ast.Add):
            return np.add(left, right, dtype=np.float32)
        if isinstance(node.op, ast.Sub):
            return np.subtract(left, right, dtype=np.float32)
        if isinstance(node.op, ast.Mult):
            return np.multiply(left, right, dtype=np.float32)
        if isinstance(node.op, ast.Div):
            return _safe_ratio(np.asarray(left), np.asarray(right))
        raise ValueError(f"不支持的二元操作: {type(node.op).__name__}")

    def visit_UnaryOp(self, node):
        v = self.visit(node.operand)
        if isinstance(node.op, ast.USub):
            return -v
        if isinstance(node.op, ast.UAdd):
            return +v
        raise ValueError(f"不支持的一元操作: {type(node.op).__name__}")

    def visit_Name(self, node):
        if node.id not in self.bn_map:
            raise ValueError(f"传感器 {self.sensor_key} 不含波段 {node.id}")
        idx = self.bn_map[node.id]
        if idx >= self.image.shape[0]:
            raise ValueError(f"波段 {node.id} 索引 {idx} 超出影像通道数 {self.image.shape[0]}")
        return self.image[idx].astype(np.float32)

    def visit_Constant(self, node):
        if not isinstance(node.value, (int, float)):
            raise ValueError("常量只允许数字")
        return float(node.value)


def eval_band_expr(image: np.ndarray, sensor: str, expr: str,
                   bn_map: Optional[Dict[str, int]] = None) -> np.ndarray:
    """
    安全求值波段表达式,如 'B5/B6'、'(B5+B7)/B6'。
    bn_map: 外部传入的 BN→通道索引映射(优先于内置常量),供 delivery_project 装载的动态数据用。
    """
    sensor_key = normalize_sensor(sensor)
    if sensor_key is None:
        raise ValueError(f"未知传感器: {sensor}")
    tree = ast.parse(expr, mode="eval")
    _validate_ast(tree)
    evaluator = _SafeBinOpEval(image, sensor_key, bn_map=bn_map)
    result = evaluator.visit(tree)
    arr = np.asarray(result, dtype=np.float32)
    if arr.ndim == 0:
        raise ValueError("表达式必须至少引用一个波段")
    return arr


# ─────────────────────────────────────────────
# Crosta PCA 法
# ─────────────────────────────────────────────

def _pc_match_score(loading: np.ndarray, band_names: List[str],
                    positive_bands: List[str], negative_bands: List[str]) -> float:
    """
    对某个 PC 的载荷向量打分:
      positive_bands 上载荷为正 → +1
      negative_bands 上载荷为负 → +1
      反向则 -1
    返回带符号得分(绝对值最大者为目标 PC,符号决定是否反号)。
    """
    name_to_load = {name: float(loading[i]) for i, name in enumerate(band_names)}
    score = 0.0
    for b in positive_bands:
        if b in name_to_load:
            score += 1.0 if name_to_load[b] > 0 else -1.0
    for b in negative_bands:
        if b in name_to_load:
            score += 1.0 if name_to_load[b] < 0 else -1.0
    return score


def calc_crosta_pca(
    image: np.ndarray,
    sensor: str,
    input_bands: List[str],
    positive_bands: List[str],
    negative_bands: List[str],
    roi_mask: Optional[np.ndarray] = None,
    bn_map: Optional[Dict[str, int]] = None,
) -> Tuple[np.ndarray, int, int, np.ndarray]:
    """
    执行 Crosta PCA。

    Returns
    -------
    index_map : (H, W) float32, ROI 外为 NaN
    target_pc_index : 0-based PC 序号(1-based 显示为 PC{idx+1})
    sign : +1 或 -1(应用到 PC score 的符号,保证"亮色=蚀变")
    loadings : (K, K) 载荷矩阵(每列是一个 PC 的载荷)
    """
    sensor_key = normalize_sensor(sensor)
    if bn_map is None:
        bn_map = _bn_map(sensor_key)

    # 取波段子集
    band_indices = []
    for b in input_bands:
        if b not in bn_map:
            raise ValueError(f"传感器 {sensor_key} 不含 PCA 输入波段 {b}")
        idx = bn_map[b]
        if idx >= image.shape[0]:
            raise ValueError(f"PCA 波段 {b} 索引 {idx} 超出影像通道数 {image.shape[0]}")
        band_indices.append(idx)

    # 组成 (H, W, K)
    subset = np.stack([image[i].astype(np.float32) for i in band_indices], axis=-1)
    H, W, K = subset.shape

    # ROI 掩膜 + 有效像素
    flat = subset.reshape(-1, K)
    finite_mask = np.all(np.isfinite(flat), axis=1)
    if roi_mask is not None:
        roi_flat = roi_mask.reshape(-1).astype(bool)
        sample_mask = finite_mask & roi_flat
    else:
        sample_mask = finite_mask

    samples = flat[sample_mask]
    if samples.shape[0] < K + 1:
        raise ValueError(f"ROI 内有效像素不足以执行 PCA (需 ≥{K+1},实际 {samples.shape[0]})")

    # 中心化 + 协方差 + 特征分解
    mean = samples.mean(axis=0)
    centered = samples - mean
    cov = np.cov(centered, rowvar=False)
    # eigh 返回升序特征值;转为降序(PC1 方差最大)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]      # 每列是一个 PC 的载荷

    # 自动选目标 PC: 匹配 pc_criterion 得分绝对值最大者
    best_pc = 0
    best_score = 0.0
    for pc in range(K):
        score = _pc_match_score(eigvecs[:, pc], input_bands, positive_bands, negative_bands)
        if abs(score) > abs(best_score):
            best_score = score
            best_pc = pc

    sign = -1 if best_score < 0 else 1

    # 把目标 PC 的 score 投影回全图
    all_centered = flat - mean
    pc_score_all = (all_centered @ eigvecs[:, best_pc]) * sign

    index_map = np.full(H * W, np.nan, dtype=np.float32)
    index_map[finite_mask] = pc_score_all[finite_mask].astype(np.float32)
    # ROI 外置 NaN(便于阈值统计也只走 ROI)
    if roi_mask is not None:
        roi_flat = roi_mask.reshape(-1).astype(bool)
        index_map[~roi_flat] = np.nan

    return index_map.reshape(H, W), best_pc, sign, eigvecs


# ─────────────────────────────────────────────
# ROI 掩膜
# ─────────────────────────────────────────────

def build_roi_mask(
    image_shape: Tuple[int, int],
    roi_geojson: Optional[Dict[str, Any]] = None,
    transform: Optional[Any] = None,
    pixel_polygon: Optional[List[List[float]]] = None,
    image_crs: Optional[Any] = None,
) -> Optional[np.ndarray]:
    """
    构建 ROI 像素掩膜 (True=ROI 内)。

    优先级:
      1. roi_geojson + transform → rasterio.features.geometry_mask
         若 image_crs 不是 EPSG:4326,先把 ROI 从 EPSG:4326 重投影到 image_crs
      2. pixel_polygon (像素坐标多边形) → PIL.ImageDraw 栅格化
      3. 都没有 → 返回 None (调用方表示"全图")
    """
    H, W = image_shape

    if roi_geojson and transform is not None:
        try:
            from rasterio.features import geometry_mask
            geom = roi_geojson["geometry"] if "geometry" in roi_geojson else roi_geojson
            # ROI 几何在 EPSG:4326;若 image CRS 不是 4326,需要 reproject
            if image_crs is not None:
                try:
                    target_epsg = image_crs.to_epsg() if hasattr(image_crs, "to_epsg") else None
                except Exception:
                    target_epsg = None
                if target_epsg and target_epsg != 4326:
                    from rasterio.warp import transform_geom
                    geom = transform_geom("EPSG:4326", image_crs, geom)
            mask = geometry_mask(
                [geom], out_shape=(H, W), transform=transform, invert=True
            )
            return mask.astype(bool)
        except Exception:
            pass

    if pixel_polygon:
        try:
            from PIL import Image as PILImage, ImageDraw
            img = PILImage.new("L", (W, H), 0)
            pts = [(float(x), float(y)) for x, y in pixel_polygon]
            ImageDraw.Draw(img).polygon(pts, outline=1, fill=1)
            return np.array(img, dtype=bool)
        except Exception:
            pass

    return None


# ─────────────────────────────────────────────
# 高光谱吸收深度 (band depth) — EnMAP 等高光谱专用
# ─────────────────────────────────────────────

def _reflectance_at(image: np.ndarray, wavelengths_um: np.ndarray,
                    target_um: float, win_um: float = 0.01) -> np.ndarray:
    """
    取 target_um 附近窗口内各波段反射率的均值(H,W)。
    窗口内无波段时退化为最近波段。窗口默认 ±10nm。
    """
    diff = np.abs(wavelengths_um - target_um)
    sel = np.where(diff <= win_um)[0]
    if sel.size == 0:
        sel = np.array([int(np.argmin(diff))])
    return np.nanmean(image[sel].astype(np.float32), axis=0)


def calc_band_depth(
    image: np.ndarray,
    wavelengths_um: np.ndarray,
    feature_um: float,
    shoulder_um: Tuple[float, float],
    sign: int = 1,
    roi_mask: Optional[np.ndarray] = None,
    win_um: float = 0.01,
) -> np.ndarray:
    """
    连续统去除吸收深度 (continuum-removed band depth)。

    在诊断吸收中心 feature_um 处,用左右肩部 shoulder_um=(left,right) 的反射率做线性连续统,
    depth = 1 - R(feature) / continuum(feature)。值越大→吸收越深→该矿物越富集。

    sign:  +1 异常=吸收深(常规,如粘土/碳酸盐/铁染);
           -1 异常=吸收浅/缺失(如红层褪色:Fe³⁺ 还原导致吸收减弱)。
    返回 index_map (H,W) float32 = sign*depth,ROI 外置 NaN。
    """
    left_um, right_um = float(shoulder_um[0]), float(shoulder_um[1])
    r_feat = _reflectance_at(image, wavelengths_um, feature_um, win_um)
    r_left = _reflectance_at(image, wavelengths_um, left_um, win_um)
    r_right = _reflectance_at(image, wavelengths_um, right_um, win_um)

    span = (right_um - left_um)
    if abs(span) < 1e-9:
        continuum = (r_left + r_right) / 2.0
    else:
        frac = (feature_um - left_um) / span
        continuum = r_left + (r_right - r_left) * frac

    with np.errstate(divide="ignore", invalid="ignore"):
        depth = 1.0 - (r_feat / continuum)
    # 连续统非正(阴影/无效)→ 无意义,置 NaN
    depth = np.where(np.isfinite(continuum) & (continuum > 1e-4), depth, np.nan)

    index_map = (float(sign) * depth).astype(np.float32)
    if roi_mask is not None:
        index_map = np.where(roi_mask, index_map, np.nan).astype(np.float32)
    return index_map


# ─────────────────────────────────────────────
# 地植物学胁迫(丛林模式)
# ─────────────────────────────────────────────

def calc_veg_stress(
    image: np.ndarray,
    sensor: str,
    roi_mask: Optional[np.ndarray] = None,
    bn_map: Optional[Dict[str, int]] = None,
    index: str = "ndre",
    veg_floor: float = 0.30,
    wavelengths_um: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    地植物学胁迫指数 —— 丛林模式的核心探测代理。

    郁闭冠层下岩石光谱被遮蔽,常规蚀变失效;但矿化/微渗漏会使其上方植物受金属毒害/
    还原胁迫 → 叶绿素↓、红边蓝移。本函数在"有植被"(NDVI>veg_floor)的像元上计算植被
    活力(默认 NDRE),再取负向 —— 活力越低输出越高,从而与 _threshold_anomaly(高值=异常)
    天然配合:输出高值 = 植被受胁迫 = 潜在矿化指示。注意此法不应被植被掩膜剔除(植被即信号)。

    支持传感器:Sentinel-2(红边 B5/B6/B7)、EnMAP/PRISMA(高光谱,需传 wavelengths_um)。
    无红边波段的传感器(Landsat/ASTER)抛 ValueError,由 analyze_batch 静默降级。

    Parameters
    ----------
    index : {"ndre","cire","rep"}  活力代理指数
    veg_floor : float  NDVI 下限,低于此值视作无冠层、不参与(置 NaN)
    wavelengths_um : 高光谱波段中心波长(µm),提供时优先按波长取红边反射率

    Returns
    -------
    index_map (H,W) float32 = -活力;非植被像元 / ROI 外置 NaN。
    """
    sensor_key = normalize_sensor(sensor) or sensor

    if wavelengths_um is not None:
        # 高光谱:按波长就近取红边附近反射率
        red = _reflectance_at(image, wavelengths_um, 0.665)
        re1 = _reflectance_at(image, wavelengths_um, 0.705)
        re2 = _reflectance_at(image, wavelengths_um, 0.740)
        re3 = _reflectance_at(image, wavelengths_um, 0.783)
        nir = _reflectance_at(image, wavelengths_um, 0.800)
    elif sensor_key == "Sentinel2":
        bm = bn_map if bn_map is not None else _bn_map("Sentinel2")

        def _b(name: str) -> np.ndarray:
            if name not in bm or bm[name] >= image.shape[0]:
                raise ValueError(f"Sentinel-2 数据缺红边波段 {name}")
            return image[bm[name]].astype(np.float32)

        red, re1, re2, re3, nir = _b("B4"), _b("B5"), _b("B6"), _b("B7"), _b("B8")
    else:
        raise ValueError(f"传感器 {sensor_key} 无红边波段,不支持地植物学胁迫法")

    ndvi = calc_ndvi(nir, red)
    if index == "rep":
        vigor = calc_rep(red, re1, re2, re3)          # nm,越大越健康(蓝移=胁迫)
    elif index == "cire":
        vigor = calc_cire(re3, re1)
    else:  # "ndre" 默认
        vigor = calc_ndre(nir, re1)

    veg = ndvi > veg_floor
    stress = (-vigor.astype(np.float32))              # 活力低 → 胁迫高
    index_map = np.where(veg & np.isfinite(stress), stress, np.nan).astype(np.float32)
    if roi_mask is not None:
        index_map = np.where(roi_mask, index_map, np.nan).astype(np.float32)
    return index_map


# ─────────────────────────────────────────────
# 阈值化
# ─────────────────────────────────────────────

def _threshold_anomaly(
    index_map: np.ndarray,
    roi_mask: Optional[np.ndarray] = None,
    method: str = "mean_std",
    k: float = 2.0,
) -> Tuple[np.ndarray, float]:
    """
    阈值化。统计样本仅取 ROI 内有效像素;输出掩膜也只在 ROI 内置 True。
    返回 (anomaly_mask, threshold)。
    """
    if roi_mask is not None:
        stat_mask = roi_mask & np.isfinite(index_map)
    else:
        stat_mask = np.isfinite(index_map)

    samples = index_map[stat_mask]
    if samples.size == 0:
        return np.zeros(index_map.shape, dtype=bool), float("nan")

    if method == "percentile":
        thr = float(np.percentile(samples, 95))
    elif method == "median_mad":
        # 稳健阈值:中位数 + k·1.4826·MAD。对拼接缝/双峰背景的抗污染能力远强于 mean+k·std。
        # 1.4826 使 MAD 在正态下等价于 σ → 同一 k 下与 mean_std 临界值近似一致。
        med = float(np.median(samples))
        mad = float(np.median(np.abs(samples - med)))
        thr = med + k * 1.4826 * mad
    else:
        thr = float(samples.mean() + k * samples.std())

    anomaly = np.isfinite(index_map) & (index_map > thr)
    if roi_mask is not None:
        anomaly &= roi_mask
    return anomaly, thr


# ─────────────────────────────────────────────
# 结果数据结构 & 主入口
# ─────────────────────────────────────────────

@dataclass
class AlterationResult:
    mineral: str
    method: str                       # "ratio" | "pca"
    sensor: str
    index_map: np.ndarray             # (H,W) float32, ROI 外 NaN
    anomaly_mask: np.ndarray          # (H,W) bool
    anomaly_ratio: float              # ROI 内异常占比 (0~1)
    threshold: float
    # method 专属信息
    ratio_expr: Optional[str] = None
    pc_used: Optional[int] = None     # 1-based PC 编号
    sign: Optional[int] = None
    warning: Optional[str] = None
    label: Optional[str] = None       # 显示名(兼容旧 API,默认 = mineral)

    def __post_init__(self):
        if self.label is None:
            self.label = self.mineral


@dataclass
class BatchResult:
    deposit_type: str
    sensor: str
    targets: List[Dict[str, Any]]                            # 每个 target 的元信息
    results: Dict[Tuple[str, str], AlterationResult]         # {(mineral, method): result}
    intersection_per_mineral: Dict[str, np.ndarray] = field(default_factory=dict)
    overall_intersection: Optional[np.ndarray] = None        # 所有矿物两方法交集的并集
    score_map: Optional[np.ndarray] = None                   # 多蚀变叠加评分
    roi_mask: Optional[np.ndarray] = None


def analyze_single(
    image: np.ndarray,
    sensor: str,
    target_spec: Dict[str, Any],
    method: str,
    roi_mask: Optional[np.ndarray] = None,
    threshold_method: str = "mean_std",
    k: float = 2.0,
    bn_map: Optional[Dict[str, int]] = None,
    wavelengths_um: Optional[np.ndarray] = None,
) -> AlterationResult:
    """
    对单个 (蚀变矿物 × 方法) 组合执行分析。

    target_spec 来自 alteration_db.get_recommended_targets():
      {mineral, ratio_expr, pca_spec, ratio_available, pca_available, ...}

    method=='band_depth' 为高光谱(EnMAP)专用,需传 wavelengths_um,
    并从 target_spec['enmap_feature'] 取 {feature_um, shoulder_um, sign}。
    """
    sensor_key = normalize_sensor(sensor) or sensor
    mineral = target_spec["mineral"]

    if method == "ratio":
        expr = target_spec.get("ratio_expr")
        if not expr or not target_spec.get("ratio_available"):
            raise ValueError(f"传感器 {sensor_key} 不支持 {mineral} 的波段比值法")
        index_map = eval_band_expr(image, sensor, expr, bn_map=bn_map)
        # ROI 外置 NaN
        if roi_mask is not None:
            index_map = np.where(roi_mask, index_map, np.nan).astype(np.float32)
        anomaly, thr = _threshold_anomaly(index_map, roi_mask, threshold_method, k)
        roi_size = int(roi_mask.sum()) if roi_mask is not None else int(np.isfinite(index_map).sum())
        anomaly_ratio = float(anomaly.sum()) / max(roi_size, 1)
        return AlterationResult(
            mineral=mineral, method="ratio", sensor=sensor_key,
            index_map=index_map, anomaly_mask=anomaly,
            anomaly_ratio=anomaly_ratio, threshold=thr,
            ratio_expr=expr,
        )

    if method == "pca":
        spec = target_spec.get("pca_spec")
        if not spec or not target_spec.get("pca_available"):
            raise ValueError(f"传感器 {sensor_key} 不支持 {mineral} 的 Crosta PCA 法")
        index_map, pc_idx, sign, _ = calc_crosta_pca(
            image, sensor,
            input_bands=spec["input_bands"],
            positive_bands=spec["positive_bands"],
            negative_bands=spec["negative_bands"],
            roi_mask=roi_mask,
            bn_map=bn_map,
        )
        anomaly, thr = _threshold_anomaly(index_map, roi_mask, threshold_method, k)
        roi_size = int(roi_mask.sum()) if roi_mask is not None else int(np.isfinite(index_map).sum())
        anomaly_ratio = float(anomaly.sum()) / max(roi_size, 1)
        return AlterationResult(
            mineral=mineral, method="pca", sensor=sensor_key,
            index_map=index_map, anomaly_mask=anomaly,
            anomaly_ratio=anomaly_ratio, threshold=thr,
            pc_used=pc_idx + 1, sign=sign,
        )

    if method == "band_depth":
        feat = target_spec.get("enmap_feature")
        if not feat or not target_spec.get("band_depth_available"):
            raise ValueError(f"传感器 {sensor_key} 不支持 {mineral} 的高光谱吸收深度法")
        if wavelengths_um is None:
            raise ValueError("band_depth 方法需要 wavelengths_um(波段中心波长)")
        feat_sign = int(feat.get("sign", 1))
        index_map = calc_band_depth(
            image, wavelengths_um,
            feature_um=float(feat["feature_um"]),
            shoulder_um=tuple(feat["shoulder_um"]),
            sign=feat_sign,
            roi_mask=roi_mask,
        )
        anomaly, thr = _threshold_anomaly(index_map, roi_mask, threshold_method, k)
        roi_size = int(roi_mask.sum()) if roi_mask is not None else int(np.isfinite(index_map).sum())
        anomaly_ratio = float(anomaly.sum()) / max(roi_size, 1)
        return AlterationResult(
            mineral=mineral, method="band_depth", sensor=sensor_key,
            index_map=index_map, anomaly_mask=anomaly,
            anomaly_ratio=anomaly_ratio, threshold=thr,
            ratio_expr=f"BD({feat['feature_um']}µm)",
            sign=feat_sign,
        )

    if method == "veg_stress":
        if not target_spec.get("veg_stress_available"):
            raise ValueError(f"传感器 {sensor_key} 不支持 {mineral} 的地植物学胁迫法")
        vspec = target_spec.get("veg_stress_spec") or {}
        index_map = calc_veg_stress(
            image, sensor,
            roi_mask=roi_mask, bn_map=bn_map,
            index=vspec.get("index", "ndre"),
            veg_floor=float(vspec.get("veg_floor", 0.30)),
            wavelengths_um=wavelengths_um,
        )
        anomaly, thr = _threshold_anomaly(index_map, roi_mask, threshold_method, k)
        roi_size = int(roi_mask.sum()) if roi_mask is not None else int(np.isfinite(index_map).sum())
        anomaly_ratio = float(anomaly.sum()) / max(roi_size, 1)
        return AlterationResult(
            mineral=mineral, method="veg_stress", sensor=sensor_key,
            index_map=index_map, anomaly_mask=anomaly,
            anomaly_ratio=anomaly_ratio, threshold=thr,
            ratio_expr=f"VEG_STRESS({vspec.get('index', 'ndre')})",
            sign=-1,
        )

    raise ValueError(f"未知方法: {method}")


def analyze_batch(
    image: np.ndarray,
    sensor: str,
    deposit_type: str,
    selected_minerals: Optional[List[str]] = None,
    methods: Optional[List[str]] = None,
    roi_mask: Optional[np.ndarray] = None,
    threshold_method: str = "mean_std",
    k: float = 2.0,
) -> BatchResult:
    """
    批量分析: 矿床类型 → 推荐蚀变矿物 × 方法。

    selected_minerals: None=用全部推荐;否则只跑列表中的矿物
    methods: None=["ratio","pca"];也可只传一个
    """
    if methods is None:
        methods = ["ratio", "pca"]

    all_targets = get_recommended_targets(deposit_type, sensor)
    if selected_minerals is not None:
        sel = set(selected_minerals)
        all_targets = [t for t in all_targets if t["mineral"] in sel]

    results: Dict[Tuple[str, str], AlterationResult] = {}
    intersection_per_mineral: Dict[str, np.ndarray] = {}
    H, W = image.shape[1], image.shape[2]
    score_map = np.zeros((H, W), dtype=np.float32)
    overall_inter = np.zeros((H, W), dtype=bool)

    for target in all_targets:
        per_method_masks: Dict[str, np.ndarray] = {}
        for m in methods:
            available = target.get(f"{m}_available", False)
            if not available:
                continue
            try:
                res = analyze_single(image, sensor, target, m, roi_mask, threshold_method, k)
                results[(target["mineral"], m)] = res
                per_method_masks[m] = res.anomaly_mask
            except Exception as e:
                results[(target["mineral"], m)] = AlterationResult(
                    mineral=target["mineral"], method=m, sensor=normalize_sensor(sensor) or sensor,
                    index_map=np.full((H, W), np.nan, dtype=np.float32),
                    anomaly_mask=np.zeros((H, W), dtype=bool),
                    anomaly_ratio=0.0, threshold=float("nan"),
                    warning=f"{m} 失败: {e}",
                )

        # 两方法交集(只在两种方法都成功时)
        if "ratio" in per_method_masks and "pca" in per_method_masks:
            inter = per_method_masks["ratio"] & per_method_masks["pca"]
            intersection_per_mineral[target["mineral"]] = inter
            overall_inter |= inter

        # 多蚀变叠加评分: priority=1 权重 1.0, priority=2 权重 0.5, priority>=3 权重 0.25
        # 分母用用户勾选的方法总数(而非实际跑成功的数),单方法天然减半 — 双方法互证才能拿满
        weight = {1: 1.0, 2: 0.5}.get(target["priority"], 0.25)
        if per_method_masks:
            stacked = sum(m.astype(np.float32) for m in per_method_masks.values())
            avg = stacked / max(len(methods), 1)
            score_map += weight * avg

    return BatchResult(
        deposit_type=deposit_type,
        sensor=normalize_sensor(sensor) or sensor,
        targets=all_targets,
        results=results,
        intersection_per_mineral=intersection_per_mineral,
        overall_intersection=overall_inter,
        score_map=score_map,
        roi_mask=roi_mask,
    )


# ─────────────────────────────────────────────
# 旧 API 兼容层 (向后兼容 /api/analyze、/api/clustering)
# ─────────────────────────────────────────────

MINERAL_CATALOG = {
    "iron_oxide": {"label": "铁氧化物", "sensors": ["Landsat8/9", "Sentinel2"],
                   "method": "ratio", "description": "Red/Blue 比值",
                   "_ratio": {"Landsat8/9": "B4/B2", "Sentinel2": "B4/B2"}},
    "aloh":       {"label": "Al-OH 羟基", "sensors": ["Landsat8/9", "Sentinel2", "ASTER"],
                   "method": "ratio", "description": "SWIR1/SWIR2",
                   "_ratio": {"Landsat8/9": "B6/B7", "Sentinel2": "B11/B12", "ASTER": "B5/B6"}},
    "carbonate":  {"label": "碳酸盐", "sensors": ["ASTER"],
                   "method": "ratio", "description": "ASTER B8/B9",
                   "_ratio": {"ASTER": "B8/B9"}},
    "ferrous":    {"label": "亚铁矿物 (绿泥石)", "sensors": ["Landsat8/9", "Sentinel2", "ASTER"],
                   "method": "ratio", "description": "NIR/Red",
                   "_ratio": {"Landsat8/9": "B5/B4", "Sentinel2": "B8/B4", "ASTER": "B4/B5"}},
    "silica":     {"label": "硅化", "sensors": ["Landsat8/9", "Sentinel2", "ASTER"],
                   "method": "ratio", "description": "SWIR1/NIR",
                   "_ratio": {"Landsat8/9": "B6/B5", "Sentinel2": "B11/B8", "ASTER": "B5/B4"}},
    "propylitic": {"label": "青磐岩化", "sensors": ["Landsat8/9", "Sentinel2", "ASTER"],
                   "method": "composite", "description": "ferrous × (1 + silica)",
                   "_ratio": {"Landsat8/9": "(B5/B4)*(1+B6/B5)",
                              "Sentinel2": "(B8/B4)*(1+B11/B8)",
                              "ASTER":     "(B4/B5)*(1+B5/B4)"}},
}


def get_supported_minerals(sensor: str) -> List[Dict]:
    """旧前端兼容: 返回固定 6 矿物的支持情况。"""
    s = sensor.replace("_L2", "")
    out = []
    for key, info in MINERAL_CATALOG.items():
        out.append({
            "key": key, "label": info["label"], "method": info["method"],
            "description": info["description"], "supported": s in info["sensors"],
        })
    return out


def analyze_alteration(
    image: np.ndarray,
    sensor: str,
    mineral: str,
    threshold_method: str = "mean_std",
    k: float = 1.5,
) -> AlterationResult:
    """旧 API 兼容: 用 MINERAL_CATALOG 中的固定表达式跑单矿物分析(无 ROI)。"""
    s = sensor.replace("_L2", "")
    if mineral not in MINERAL_CATALOG:
        raise ValueError(f"未知矿物类型: {mineral}")
    info = MINERAL_CATALOG[mineral]
    if s not in info["_ratio"]:
        raise ValueError(f"传感器 {s} 不支持矿物 {info['label']}")

    target = {
        "mineral": info["label"],
        "ratio_expr": info["_ratio"][s],
        "ratio_available": True,
        "pca_available": False,
    }
    return analyze_single(image, s, target, "ratio", None, threshold_method, k)


if __name__ == "__main__":
    # 演示: 用合成 Landsat 数据
    rng = np.random.default_rng(42)
    image = rng.uniform(0.05, 0.3, size=(7, 64, 64)).astype(np.float32)
    image[5, 24:40, 24:40] *= 2.0  # SWIR1 异常

    print("=== 新 API: analyze_batch ===")
    batch = analyze_batch(image, "Landsat8/9", "斑岩型铜矿")
    for (mineral, method), res in batch.results.items():
        info = f"{res.threshold:.3f}" if np.isfinite(res.threshold) else "n/a"
        extra = ""
        if method == "pca" and res.pc_used:
            extra = f" PC{res.pc_used} sign={'+' if res.sign>0 else '-'}"
        warn = f"  ⚠ {res.warning}" if res.warning else ""
        print(f"  {mineral:<10} [{method:<5}] 异常占比 {res.anomaly_ratio*100:5.2f}% 阈值={info}{extra}{warn}")
    print(f"  综合交集像素: {int(batch.overall_intersection.sum())}")
    print(f"  评分图 max={batch.score_map.max():.2f}")

    print()
    print("=== 旧 API: analyze_alteration ===")
    res = analyze_alteration(image, "Landsat8/9_L2", "aloh")
    print(f"  {res.mineral} 异常占比 {res.anomaly_ratio*100:.2f}%")
