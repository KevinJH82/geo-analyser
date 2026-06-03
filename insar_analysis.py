"""
insar_analysis.py — geo-analyser InSAR 异常检测模块(Phase 1.5)

提供三个核心函数:
- coherence_to_stability(): 相干性 → 地表稳定性得分
- los_velocity_clustering(): LOS 速率聚类识别活跃形变区
- fusion_deformation_mineral(): 形变-矿物联合异常分析

与现有 alteration_analysis.py 的 API 风格保持一致(输入 ndarray + 阈值,返回结果 dataclass)。
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


# 相干性的几个语义阈值(经验值,可以在 UI 配置)
COHERENCE_HIGH = 0.6      # 高度可靠,稳定区域
COHERENCE_MID = 0.3       # 中等,可用
COHERENCE_LOW = 0.15      # 低,接近噪声


@dataclass
class InsarAnalysisResult:
    """InSAR 分析结果(与 alteration_analysis.AlterationResult 风格一致)。"""
    name: str
    array: np.ndarray                     # 主结果数组(float32)
    mask: Optional[np.ndarray] = None     # 二值掩膜(True 表示异常/活跃)
    stats: Dict = field(default_factory=dict)
    colormap: str = "RdBu_r"              # 默认双向 diverging(形变方向)
    vmin: Optional[float] = None
    vmax: Optional[float] = None
    unit: str = ""
    error: Optional[str] = None


def coherence_to_stability(coherence: np.ndarray,
                           threshold: float = COHERENCE_MID) -> InsarAnalysisResult:
    """
    把相干性图转成"地表稳定性"得分(0-1),并标记滑坡易发区(低相干 + 非水体)。

    语义:
    - 高相干(>0.6):稳定 → 高分
    - 低相干(<0.15):可能是植被、水体、剧烈形变 → 低分
    - 中间值:线性映射

    Parameters
    ----------
    coherence : np.ndarray, 0-1
    threshold : float, 低于此值标记为"不稳定/异常",默认 0.3

    Returns
    -------
    InsarAnalysisResult
    """
    coh = np.asarray(coherence, dtype=np.float32)
    coh_clipped = np.clip(coh, 0, 1)

    # 简单线性映射:相干性即稳定性得分
    stability = coh_clipped

    # 异常掩膜:相干性 < threshold
    mask = (coh_clipped < threshold) & np.isfinite(coh_clipped)

    stats = {
        "mean": float(np.nanmean(coh_clipped)),
        "median": float(np.nanmedian(coh_clipped)),
        "low_coherence_ratio": float(np.sum(mask) / mask.size) if mask.size else 0.0,
        "threshold": threshold,
    }
    return InsarAnalysisResult(
        name="coherence_stability",
        array=stability,
        mask=mask,
        stats=stats,
        colormap="RdYlGn",
        vmin=0.0, vmax=1.0,
        unit="coherence (0-1)",
    )


def los_velocity_clustering(velocity: np.ndarray,
                            coherence: Optional[np.ndarray] = None,
                            velocity_threshold_mm_year: float = 5.0,
                            coherence_threshold: float = COHERENCE_MID,
                            eps_pixels: int = 5,
                            min_samples: int = 50) -> InsarAnalysisResult:
    """
    LOS 速率聚类:识别形变活跃区(超过阈值的连通块)。

    使用 DBSCAN-like 简化算法(连通块标记 + 大小过滤),无需 scikit-learn。

    Parameters
    ----------
    velocity : np.ndarray, mm/year
    coherence : np.ndarray, optional, 用于掩膜
    velocity_threshold_mm_year : 形变量阈值,绝对值超过即视为活跃像素
    coherence_threshold : 相干性阈值
    eps_pixels : 形态学闭运算的核大小(填洞用)
    min_samples : 最小连通块像素数

    Returns
    -------
    InsarAnalysisResult, .array 是活跃区掩膜(0/1 float32),.stats 含聚类数等
    """
    from scipy.ndimage import binary_closing, label

    v = np.asarray(velocity, dtype=np.float32)
    if coherence is not None:
        mask_v = (np.abs(v) >= velocity_threshold_mm_year) & (np.asarray(coherence) >= coherence_threshold)
    else:
        mask_v = np.abs(v) >= velocity_threshold_mm_year

    # 形态学闭运算填洞
    mask_closed = binary_closing(mask_v, structure=np.ones((eps_pixels, eps_pixels)))

    # 标记连通块,过滤太小的
    labeled, n_labels = label(mask_closed.astype(np.uint8))
    final = np.zeros_like(mask_closed, dtype=np.uint8)
    keep = 0
    for i in range(1, n_labels + 1):
        size = int(np.sum(labeled == i))
        if size >= min_samples:
            final[labeled == i] = 1
            keep += 1

    # 区分沉降区(负速率)/抬升区(正速率)
    subside_ratio = float(np.sum((final == 1) & (v < 0)) / max(1, int(np.sum(final == 1))))

    stats = {
        "active_cluster_count": keep,
        "active_pixel_ratio": float(final.mean()) if final.size else 0.0,
        "subsidence_ratio_within_active": subside_ratio,
        "velocity_threshold_mm_year": velocity_threshold_mm_year,
        "coherence_threshold": coherence_threshold,
    }
    return InsarAnalysisResult(
        name="los_velocity_clusters",
        array=final.astype(np.float32),
        mask=(final == 1),
        stats=stats,
        colormap="hot",
        vmin=0.0, vmax=1.0,
        unit="active region mask",
    )


def fusion_deformation_mineral(velocity: np.ndarray,
                               mineral_anomaly: np.ndarray,
                               coherence: Optional[np.ndarray] = None,
                               weight_v: float = 0.5,
                               weight_m: float = 0.5,
                               coherence_threshold: float = COHERENCE_MID) -> InsarAnalysisResult:
    """
    形变-矿物联合异常分析。

    场景:在矿区周边,如果同一位置既有 LOS 速率异常(活跃形变),又有蚀变矿物异常,
    则该位置很可能是活跃开采区或活动构造控矿区。

    Parameters
    ----------
    velocity : np.ndarray, LOS 速率(mm/year)
    mineral_anomaly : np.ndarray, 矿物异常得分(任意尺度,会标准化)
    coherence : np.ndarray, optional
    weight_v, weight_m : 形变权重 vs 矿物权重(应和为 1.0)
    coherence_threshold : 相干性阈值

    Returns
    -------
    InsarAnalysisResult, .array 是联合异常得分(0-1)
    """
    v = np.abs(np.asarray(velocity, dtype=np.float32))
    m = np.asarray(mineral_anomaly, dtype=np.float32)

    # 标准化到 [0,1] 用 percentile 截断防极端值
    def _normalize(arr, p_lo=2, p_hi=98):
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return np.zeros_like(arr)
        lo = np.percentile(finite, p_lo)
        hi = np.percentile(finite, p_hi)
        if hi - lo < 1e-9:
            return np.zeros_like(arr)
        out = np.clip((arr - lo) / (hi - lo), 0, 1)
        return out

    v_n = _normalize(v)
    m_n = _normalize(m)
    fused = weight_v * v_n + weight_m * m_n

    # 应用相干性掩膜
    if coherence is not None:
        mask_coh = np.asarray(coherence) >= coherence_threshold
        fused = np.where(mask_coh, fused, 0.0)

    # 高分阈值(top 5%)作为强联合异常
    finite = fused[np.isfinite(fused)]
    high_thr = float(np.percentile(finite, 95)) if finite.size else 1.0
    mask = fused >= high_thr

    stats = {
        "high_anomaly_threshold": high_thr,
        "high_anomaly_pixel_ratio": float(mask.mean()) if mask.size else 0.0,
        "weight_velocity": weight_v,
        "weight_mineral": weight_m,
    }
    return InsarAnalysisResult(
        name="deformation_mineral_fusion",
        array=fused,
        mask=mask,
        stats=stats,
        colormap="hot",
        vmin=0.0, vmax=1.0,
        unit="联合异常得分 (0-1)",
    )
