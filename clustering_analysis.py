"""
尖点突破分析模块
基于空间聚类算法识别蚀变异常集中形成的突破区
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# ─────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────

@dataclass
class ClusterStats:
    """单个聚类的统计信息"""
    cluster_id: int
    pixel_count: int          # 像元数量
    density_score: float      # 局部密度得分 [0,1]
    mean_index: float         # 聚类内指数均值
    max_index: float          # 聚类内指数最大值
    centroid: Tuple[int, int] # 质心 (row, col)
    bbox: Tuple[int, int, int, int]  # (min_row, min_col, max_row, max_col)


@dataclass
class BreakthroughZone:
    """突破区"""
    cluster_id: int
    pixel_count: int
    density_score: float      # 密度得分（越高越集中）
    breakthrough_strength: float  # 突破强度 = 密度得分 × 均值指数
    centroid: Tuple[int, int]
    bbox: Tuple[int, int, int, int]
    mean_index: float
    max_index: float


@dataclass
class ClusteringResult:
    """聚类分析完整结果"""
    mineral: str
    label: str
    sensor: str
    algorithm: str             # 使用的算法
    n_clusters: int            # 有效聚类数
    cluster_label_map: np.ndarray   # 标签图（-1=噪声，0,1,2,...=聚类ID）
    density_map: np.ndarray    # 密度热力图（float32，归一化到 [0,1]）
    cluster_stats: List[ClusterStats]
    breakthrough_zones: List[BreakthroughZone]  # 突破区列表（按强度降序）
    spatial_autocorrelation: float   # Moran's I 空间自相关
    warning: Optional[str] = None


# ─────────────────────────────────────────────
# 密度计算
# ─────────────────────────────────────────────

def _gaussian_density(
    anomaly_mask: np.ndarray,
    sigma: float = 5.0,
) -> np.ndarray:
    """
    用高斯核计算局部异常密度图。
    返回归一化到 [0,1] 的 float32 数组。
    """
    from scipy.ndimage import gaussian_filter
    density = gaussian_filter(anomaly_mask.astype(np.float32), sigma=sigma)
    dmax = density.max()
    if dmax > 0:
        density /= dmax
    return density


# ─────────────────────────────────────────────
# DBSCAN 空间聚类
# ─────────────────────────────────────────────

def _dbscan_cluster(
    anomaly_mask: np.ndarray,
    eps: float = 10.0,
    min_samples: int = 5,
) -> np.ndarray:
    """
    对异常像元做 DBSCAN 空间聚类。

    Parameters
    ----------
    anomaly_mask : bool ndarray (H, W)
    eps : float   邻域半径（像元单位）
    min_samples : int  核心点最小邻居数

    Returns
    -------
    label_map : int32 ndarray (H, W)
        -1=噪声，0,1,... = 聚类 ID
    """
    from sklearn.cluster import DBSCAN

    # 提取异常像元坐标
    rows, cols = np.where(anomaly_mask)
    if len(rows) < min_samples:
        return np.full(anomaly_mask.shape, -1, dtype=np.int32)

    coords = np.column_stack([rows, cols]).astype(np.float32)
    labels = DBSCAN(eps=eps, min_samples=min_samples, algorithm='ball_tree',
                    metric='euclidean', n_jobs=-1).fit_predict(coords)

    label_map = np.full(anomaly_mask.shape, -1, dtype=np.int32)
    label_map[rows, cols] = labels
    return label_map


# ─────────────────────────────────────────────
# 层次聚类
# ─────────────────────────────────────────────

def _hierarchical_cluster(
    anomaly_mask: np.ndarray,
    n_clusters: int = 5,
    max_pixels: int = 5000,
) -> np.ndarray:
    """
    对异常像元做层次聚类（Ward linkage）。
    像元过多时随机采样后反映射。
    """
    from sklearn.cluster import AgglomerativeClustering
    from scipy.spatial import KDTree

    rows, cols = np.where(anomaly_mask)
    if len(rows) < n_clusters:
        return np.full(anomaly_mask.shape, -1, dtype=np.int32)

    coords = np.column_stack([rows, cols]).astype(np.float32)

    # 像元过多则采样
    if len(coords) > max_pixels:
        idx = np.random.choice(len(coords), max_pixels, replace=False)
        sample_coords = coords[idx]
        sample_labels = AgglomerativeClustering(
            n_clusters=n_clusters, linkage='ward'
        ).fit_predict(sample_coords)
        # 用 KDTree 把标签映射回全部像元
        tree = KDTree(sample_coords)
        _, nearest = tree.query(coords, k=1)
        labels = sample_labels[nearest]
    else:
        labels = AgglomerativeClustering(
            n_clusters=n_clusters, linkage='ward'
        ).fit_predict(coords)

    label_map = np.full(anomaly_mask.shape, -1, dtype=np.int32)
    label_map[rows, cols] = labels
    return label_map


# ─────────────────────────────────────────────
# 聚类统计
# ─────────────────────────────────────────────

def _compute_cluster_stats(
    label_map: np.ndarray,
    index_map: np.ndarray,
    density_map: np.ndarray,
) -> List[ClusterStats]:
    stats = []
    unique_ids = [i for i in np.unique(label_map) if i >= 0]
    for cid in unique_ids:
        mask = label_map == cid
        rr, cc = np.where(mask)
        vals = index_map[mask]
        valid_vals = vals[np.isfinite(vals)]
        density_vals = density_map[mask]

        centroid = (int(rr.mean()), int(cc.mean()))
        bbox = (int(rr.min()), int(cc.min()), int(rr.max()), int(cc.max()))

        stats.append(ClusterStats(
            cluster_id=int(cid),
            pixel_count=int(mask.sum()),
            density_score=float(density_vals.mean()),
            mean_index=float(valid_vals.mean()) if valid_vals.size > 0 else 0.0,
            max_index=float(valid_vals.max()) if valid_vals.size > 0 else 0.0,
            centroid=centroid,
            bbox=bbox,
        ))
    return stats


# ─────────────────────────────────────────────
# 突破区识别
# ─────────────────────────────────────────────

def _identify_breakthrough_zones(
    cluster_stats: List[ClusterStats],
    density_threshold: float = 0.3,
    min_pixels: int = 10,
) -> List[BreakthroughZone]:
    """
    从聚类统计中筛选并量化突破区。
    突破强度 = density_score × mean_index（越高越显著）
    """
    zones = []
    for cs in cluster_stats:
        if cs.density_score < density_threshold or cs.pixel_count < min_pixels:
            continue
        zones.append(BreakthroughZone(
            cluster_id=cs.cluster_id,
            pixel_count=cs.pixel_count,
            density_score=cs.density_score,
            breakthrough_strength=cs.density_score * cs.mean_index,
            centroid=cs.centroid,
            bbox=cs.bbox,
            mean_index=cs.mean_index,
            max_index=cs.max_index,
        ))
    # 按突破强度降序
    zones.sort(key=lambda z: z.breakthrough_strength, reverse=True)
    return zones


# ─────────────────────────────────────────────
# Moran's I 空间自相关
# ─────────────────────────────────────────────

def _morans_i(anomaly_mask: np.ndarray) -> float:
    """
    计算 Moran's I 空间自相关指数（基于 4-邻域权重）。
    结果范围 [-1, 1]，>0 表示正相关（聚集），<0 表示分散。
    """
    x = anomaly_mask.astype(np.float32)
    n = x.size
    x_mean = x.mean()
    x_dev = x - x_mean

    # 4-邻域差积
    lag = (
        np.roll(x_dev, 1, axis=0) +
        np.roll(x_dev, -1, axis=0) +
        np.roll(x_dev, 1, axis=1) +
        np.roll(x_dev, -1, axis=1)
    )
    W = 4 * n  # 总权重
    numerator = n * float(np.sum(x_dev * lag))
    denominator = W * float(np.sum(x_dev ** 2))
    if denominator == 0:
        return 0.0
    return float(numerator / denominator)


# ─────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────

def analyze_clustering(
    anomaly_mask: np.ndarray,
    index_map: np.ndarray,
    mineral: str,
    label: str,
    sensor: str,
    algorithm: str = "dbscan",
    # DBSCAN 参数
    eps: float = 10.0,
    min_samples: int = 5,
    # 层次聚类参数
    n_clusters: int = 5,
    # 密度参数
    density_sigma: float = 5.0,
    # 突破区阈值
    density_threshold: float = 0.3,
    min_pixels: int = 10,
) -> ClusteringResult:
    """
    对蚀变异常掩膜执行空间聚类分析，识别突破区。

    Parameters
    ----------
    anomaly_mask : bool ndarray (H, W)
    index_map    : float32 ndarray (H, W)  连续指数值
    mineral      : 矿物 key
    label        : 显示名称
    sensor       : 传感器名称
    algorithm    : "dbscan" 或 "hierarchical"
    """
    # 聚类
    if algorithm == "hierarchical":
        label_map = _hierarchical_cluster(anomaly_mask, n_clusters=n_clusters)
    else:
        label_map = _dbscan_cluster(anomaly_mask, eps=eps, min_samples=min_samples)

    unique_ids = [i for i in np.unique(label_map) if i >= 0]
    n_clusters_found = len(unique_ids)

    # 密度图
    density_map = _gaussian_density(anomaly_mask, sigma=density_sigma)

    # 聚类统计
    cluster_stats = _compute_cluster_stats(label_map, index_map, density_map)

    # 突破区
    breakthrough_zones = _identify_breakthrough_zones(
        cluster_stats, density_threshold=density_threshold, min_pixels=min_pixels
    )

    # 空间自相关
    moran = _morans_i(anomaly_mask)

    warning = None
    if n_clusters_found == 0:
        warning = "未检测到有效聚类，请降低最小样本数或检查异常掩膜"

    return ClusteringResult(
        mineral=mineral,
        label=label,
        sensor=sensor,
        algorithm=algorithm,
        n_clusters=n_clusters_found,
        cluster_label_map=label_map,
        density_map=density_map,
        cluster_stats=cluster_stats,
        breakthrough_zones=breakthrough_zones,
        spatial_autocorrelation=moran,
        warning=warning,
    )
