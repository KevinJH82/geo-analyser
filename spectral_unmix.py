"""
spectral_unmix.py — 亚像元光谱解混 + 天然出露靶向(丛林模式增强)

丛林区两个抢救思路:
  1. 解混(unmix):半植被像元的反射率是 {绿色植被, 土壤/岩石, 阴影} 的混合。用 NNLS 线性
     解混恢复各端元丰度,把"岩石分量丰度"较高的像元拣选出来,而不是被植被掩膜一刀切丢弃。
  2. 天然出露靶向:用 BSI(裸土指数)+ 低 NDVI 定位河岸/山脊/滑坡疤/采坑等少数裸露像元,
     在这些像元上做的常规蚀变才可信 —— 据此对蚀变置信度加权。

两者均为只读纯函数,带合成数据单测,可被 app 的丛林模式按需调用。
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# 纯指数函数从 commons 共享层引入(geo-analyser 不依赖 geo-preprocess)
import os as _os, sys as _sys
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _REPO_ROOT not in _sys.path:
    _sys.path.insert(0, _REPO_ROOT)
from commons.spectral_indices import calc_ndvi, calc_bsi, _safe_divide

# Landsat8/9 波段顺序索引(0-based): B1..B7 = Coastal,Blue,Green,Red,NIR,SWIR1,SWIR2
# 仅供解混自动端元估计的回退路径与单测使用(调用方通常显式传波段)
_LANDSAT_IDX = {"blue": 1, "green": 2, "red": 3, "nir": 4, "swir1": 5, "swir2": 6}


# ─────────────────────────────────────────────
# 线性光谱解混(NNLS)
# ─────────────────────────────────────────────

@dataclass
class UnmixResult:
    abundances: np.ndarray              # (n_endmembers, H, W) float32,各端元丰度 [0,1]
    endmember_names: List[str]
    rock_fraction: np.ndarray           # (H,W) 土壤/岩石端元丰度(便捷别名)
    rmse: np.ndarray                    # (H,W) 重建残差(拟合质量)
    stats: Dict = field(default_factory=dict)


def estimate_endmembers(
    image: np.ndarray,
    ndvi: np.ndarray,
    bsi: np.ndarray,
    roi_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, List[str]]:
    """
    从影像极值自动估计三端元光谱(无需外部光谱库):
      - 绿色植被:NDVI 最高的一批像元均值
      - 土壤/岩石:BSI 最高(裸土)像元均值
      - 阴影:总反射率最低像元均值

    返回 (E (n_bands, 3), names)。
    """
    B = image.shape[0]
    flat = image.reshape(B, -1).astype(np.float32)
    nd = ndvi.ravel()
    bs = bsi.ravel()
    brightness = np.nanmean(flat, axis=0)

    valid = np.isfinite(nd) & np.isfinite(bs) & np.isfinite(brightness)
    if roi_mask is not None:
        valid &= roi_mask.ravel()
    if valid.sum() < 30:
        valid = np.isfinite(brightness)

    def _top_mean(score, frac=0.02, invert=False):
        s = score.copy()
        s[~valid] = np.nan if not invert else np.nan
        order = np.argsort(s)
        order = order[np.isfinite(s[order])]
        if order.size == 0:
            return np.nanmean(flat[:, valid], axis=1)
        n = max(5, int(order.size * frac))
        pick = order[-n:] if not invert else order[:n]
        return np.nanmean(flat[:, pick], axis=1)

    e_veg  = _top_mean(nd)                       # NDVI 高 → 植被
    e_rock = _top_mean(bs)                       # BSI 高 → 裸土/岩石
    e_shade = _top_mean(brightness, invert=True) # 亮度低 → 阴影
    E = np.stack([e_veg, e_rock, e_shade], axis=1).astype(np.float32)  # (B,3)
    return E, ["vegetation", "rock_soil", "shade"]


def linear_unmix(
    image: np.ndarray,
    endmembers: Optional[np.ndarray] = None,
    endmember_names: Optional[List[str]] = None,
    roi_mask: Optional[np.ndarray] = None,
    candidate_mask: Optional[np.ndarray] = None,
    sum_to_one: bool = True,
) -> UnmixResult:
    """
    线性光谱解混: image(b) ≈ Σ_k a_k · E(b,k), a_k ≥ 0。用 scipy NNLS 逐像元求解。

    为控制耗时,只在 candidate_mask(默认=ROI 内)像元上解混,其余置 NaN。
    endmembers 为 None 时用 estimate_endmembers 自动估计(需可计算 NDVI/BSI 的多光谱影像,
    其波段顺序假定为 Landsat8/9(见 _LANDSAT_IDX),由调用方保证)。

    sum_to_one: True 时对 NNLS 解做和归一(丰度比例),便于解释。
    """
    from scipy.optimize import nnls

    B, H, W = image.shape
    if endmembers is None:
        # 需要 NDVI/BSI;调用方应传 6 波段(含 nir/swir1)的多光谱。这里按 Landsat 顺序兜底。
        _LB = _LANDSAT_IDX
        nir = image[_LB["nir"]]; red = image[_LB["red"]]
        blue = image[_LB["blue"]]; swir1 = image[_LB["swir1"]]
        ndvi = calc_ndvi(nir, red)
        bsi = calc_bsi(blue, red, nir, swir1)
        endmembers, endmember_names = estimate_endmembers(image, ndvi, bsi, roi_mask)
    if endmember_names is None:
        endmember_names = [f"em{i}" for i in range(endmembers.shape[1])]

    K = endmembers.shape[1]
    abundances = np.full((K, H, W), np.nan, dtype=np.float32)
    rmse = np.full((H, W), np.nan, dtype=np.float32)

    if candidate_mask is None:
        candidate_mask = roi_mask if roi_mask is not None else np.ones((H, W), dtype=bool)
    cand = candidate_mask & np.all(np.isfinite(image), axis=0)

    idx = np.argwhere(cand)
    E = endmembers.astype(np.float64)
    for (r, c) in idx:
        y = image[:, r, c].astype(np.float64)
        try:
            a, res = nnls(E, y)
        except Exception:
            continue
        if sum_to_one and a.sum() > 1e-9:
            a = a / a.sum()
        abundances[:, r, c] = a.astype(np.float32)
        recon = E @ a
        rmse[r, c] = float(np.sqrt(np.mean((recon - y) ** 2)))

    rock_idx = endmember_names.index("rock_soil") if "rock_soil" in endmember_names else 1
    rock_fraction = abundances[rock_idx]

    stats = {
        "n_unmixed": int(cand.sum()),
        "mean_rmse": float(np.nanmean(rmse)) if np.isfinite(rmse).any() else None,
        "endmembers": endmember_names,
    }
    return UnmixResult(
        abundances=abundances, endmember_names=endmember_names,
        rock_fraction=rock_fraction, rmse=rmse, stats=stats,
    )


# ─────────────────────────────────────────────
# 天然出露靶向
# ─────────────────────────────────────────────

@dataclass
class ExposureResult:
    mask: np.ndarray                    # (H,W) bool,True=天然出露(裸岩/裸土)
    score: np.ndarray                   # (H,W) [0,1] 出露程度连续分
    stats: Dict = field(default_factory=dict)


def detect_natural_exposures(
    nir: np.ndarray, red: np.ndarray,
    blue: np.ndarray, swir1: np.ndarray,
    roi_mask: Optional[np.ndarray] = None,
    ndvi_max: float = 0.25,
    bsi_min: float = 0.0,
) -> ExposureResult:
    """
    检测天然出露像元(河岸/山脊/滑坡疤/路堑/采坑等少数裸露):
      出露 = (NDVI < ndvi_max) 且 (BSI > bsi_min)。
    score = 归一化 BSI · (1 - 归一化 NDVI),连续刻画"裸露程度",供蚀变置信度加权。
    丛林区出露像元少而珍贵,在其上做的常规蚀变才可信。
    """
    ndvi = calc_ndvi(nir, red)
    bsi = calc_bsi(blue, red, nir, swir1)

    mask = (ndvi < ndvi_max) & (bsi > bsi_min) & np.isfinite(ndvi) & np.isfinite(bsi)
    if roi_mask is not None:
        mask &= roi_mask

    def _norm01(a):
        f = a[np.isfinite(a)]
        if f.size == 0:
            return np.zeros_like(a, dtype=np.float32)
        lo, hi = np.percentile(f, 2), np.percentile(f, 98)
        if hi - lo < 1e-9:
            return np.zeros_like(a, dtype=np.float32)
        return np.clip((a - lo) / (hi - lo), 0, 1).astype(np.float32)

    score = (_norm01(bsi) * (1.0 - _norm01(ndvi))).astype(np.float32)
    if roi_mask is not None:
        score = np.where(roi_mask, score, np.nan).astype(np.float32)

    denom = int(roi_mask.sum()) if roi_mask is not None else mask.size
    stats = {
        "exposure_pixels": int(mask.sum()),
        "exposure_ratio": float(mask.sum()) / max(denom, 1),
        "ndvi_max": ndvi_max, "bsi_min": bsi_min,
    }
    return ExposureResult(mask=mask, score=score, stats=stats)


# ─────────────────────────────────────────────
# 演示 / 单元测试
# ─────────────────────────────────────────────

def _demo():
    np.random.seed(11)
    # 7 波段(Landsat 顺序: B1..B7 = Coastal,Blue,Green,Red,NIR,SWIR1,SWIR2)
    LB = _LANDSAT_IDX
    H, W = 50, 50
    img = np.random.uniform(0.02, 0.08, (7, H, W)).astype(np.float32)
    # 植被区(左半): 高 NIR 低 red
    img[LB["nir"], :, :25] = 0.50; img[LB["red"], :, :25] = 0.04
    # 裸岩出露(右半): 高 SWIR/red, 低 NIR
    img[LB["swir1"], :, 25:] = 0.40; img[LB["red"], :, 25:] = 0.30; img[LB["nir"], :, 25:] = 0.25
    roi = np.ones((H, W), dtype=bool)

    # 出露检测
    er = detect_natural_exposures(img[LB["nir"]], img[LB["red"]], img[LB["blue"]], img[LB["swir1"]], roi_mask=roi)
    left = er.mask[:, :25].mean(); right = er.mask[:, 25:].mean()
    print("出露检测  植被半区命中率: %.3f   裸岩半区命中率: %.3f" % (left, right))
    assert right > 0.8 and left < 0.05, "应只在裸岩半区判为出露"

    # 解混: 裸岩半区 rock_fraction 应高
    ur = linear_unmix(img, roi_mask=roi)
    rf_left = np.nanmean(ur.rock_fraction[:, :25])
    rf_right = np.nanmean(ur.rock_fraction[:, 25:])
    print("解混  植被半区岩石丰度: %.3f   裸岩半区岩石丰度: %.3f" % (rf_left, rf_right))
    print("端元:", ur.endmember_names, "| 平均 RMSE: %.4f" % (ur.stats["mean_rmse"] or 0))
    assert rf_right > rf_left, "裸岩半区岩石端元丰度应更高"
    print("PASS: spectral_unmix 工作正常")


if __name__ == "__main__":
    _demo()
