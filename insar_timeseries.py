"""
insar_timeseries.py — geo-analyser InSAR 时序分析模块(Phase 1.5)

提供:
- load_insar_stack(): 加载 geo-insar 标准输出目录下的多对 InSAR
- temporal_velocity_trend(): 多对形变量 → 线性速率(类 PS/SBAS 简化版)
- coherence_decay_model(): 相干性随时间基线的衰减建模

不依赖 MintPy/StaMPS,只用 numpy/rasterio。
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


def _read_geotiff(path: str) -> Tuple[np.ndarray, Dict]:
    """读 GeoTIFF,返回 (array, meta)。"""
    import rasterio
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        meta = {
            "shape": arr.shape,
            "transform": src.transform,
            "crs": str(src.crs) if src.crs else None,
            "nodata": src.nodata,
        }
    return arr, meta


def load_insar_stack(aoi_output_dir: str) -> Dict:
    """
    扫描 geo-insar 输出的 AOI 目录,加载所有干涉对的 LOS 形变 + 相干性。

    目录契约(commons/insar_schema.json):
      {aoi}/sentinel1_insar/<refdate>_<secdate>_<pol>/
        ├── los_displacement.tif
        ├── coherence.tif
        └── metadata.json

    Returns
    -------
    {
        "aoi_name": str,
        "n_pairs": int,
        "ref_shape": (H, W),
        "displacements": [(ref_date, sec_date, disp_array), ...],
        "coherences": [coh_array, ...],
        "metas": [pair_metadata, ...]
    }
    """
    aoi_path = Path(aoi_output_dir)
    if not aoi_path.exists():
        raise FileNotFoundError(f"AOI 目录不存在: {aoi_output_dir}")

    pair_dirs = sorted(aoi_path.glob("sentinel1_insar/*"))
    pair_dirs = [p for p in pair_dirs if p.is_dir() and (p / "metadata.json").exists()]
    if not pair_dirs:
        return {"aoi_name": aoi_path.name, "n_pairs": 0,
                "ref_shape": None, "displacements": [], "coherences": [], "metas": []}

    displacements, coherences, metas = [], [], []
    ref_shape = None

    for pdir in pair_dirs:
        with open(pdir / "metadata.json", "r", encoding="utf-8") as f:
            meta = json.load(f)
        metas.append(meta)
        disp_file = pdir / "los_displacement.tif"
        coh_file = pdir / "coherence.tif"
        if not disp_file.exists():
            continue
        disp, m = _read_geotiff(str(disp_file))
        ref_shape = ref_shape or disp.shape
        # 简单 align(如果有 shape 不一致,用 nearest 重采样)
        if disp.shape != ref_shape:
            from scipy.ndimage import zoom
            zy = ref_shape[0] / disp.shape[0]
            zx = ref_shape[1] / disp.shape[1]
            disp = zoom(disp, (zy, zx), order=1)
        displacements.append((meta["master_date"], meta["slave_date"], disp))

        if coh_file.exists():
            coh, _ = _read_geotiff(str(coh_file))
            if coh.shape != ref_shape:
                from scipy.ndimage import zoom
                zy = ref_shape[0] / coh.shape[0]
                zx = ref_shape[1] / coh.shape[1]
                coh = zoom(coh, (zy, zx), order=1)
            coherences.append(coh)
        else:
            coherences.append(None)

    return {
        "aoi_name": aoi_path.name,
        "n_pairs": len(displacements),
        "ref_shape": ref_shape,
        "displacements": displacements,
        "coherences": coherences,
        "metas": metas,
    }


def temporal_velocity_trend(stack: Dict, min_coherence: float = 0.3) -> Dict:
    """
    多对形变量 → 线性速率(mm/year)。

    简化的 SBAS 思路:
    1. 把每对 (date1, date2, disp) 转换成 (t_mid, disp / dt_years)
    2. 像素级取多对的中位数速率,鲁棒于个别对的噪声

    比真正 PS/SBAS 简陋很多,但对小 AOI 演示足够。

    Parameters
    ----------
    stack : load_insar_stack() 返回
    min_coherence : 像素级筛选阈值

    Returns
    -------
    {
        "velocity_mm_per_year": ndarray,  # 速率图
        "n_valid_pairs_per_pixel": ndarray,  # 每个像素参与中位数计算的对数
        "stats": {"mean_velocity", "median_velocity", "subsidence_max", "uplift_max"}
    }
    """
    if stack["n_pairs"] == 0:
        return {"velocity_mm_per_year": None, "n_valid_pairs_per_pixel": None, "stats": {}}

    ref_shape = stack["ref_shape"]
    pair_rates = []  # 每对一个 (rate_array, weight_mask)
    for (d1, d2, disp), coh in zip(stack["displacements"], stack["coherences"]):
        t1 = datetime.strptime(d1, "%Y-%m-%d")
        t2 = datetime.strptime(d2, "%Y-%m-%d")
        dt_days = (t2 - t1).days
        if dt_days <= 0:
            continue
        rate = disp * (365.25 / dt_days)  # mm/year
        valid = np.isfinite(rate)
        if coh is not None:
            valid &= coh >= min_coherence
        rate_masked = np.where(valid, rate, np.nan)
        pair_rates.append(rate_masked)

    if not pair_rates:
        return {"velocity_mm_per_year": None, "n_valid_pairs_per_pixel": None, "stats": {}}

    rates_3d = np.stack(pair_rates, axis=0)  # (N, H, W)
    median_rate = np.nanmedian(rates_3d, axis=0)
    n_valid = np.sum(np.isfinite(rates_3d), axis=0)

    finite = median_rate[np.isfinite(median_rate)]
    stats = {
        "mean_velocity": float(np.mean(finite)) if finite.size else 0.0,
        "median_velocity": float(np.median(finite)) if finite.size else 0.0,
        "subsidence_max_mm_year": float(np.min(finite)) if finite.size else 0.0,
        "uplift_max_mm_year": float(np.max(finite)) if finite.size else 0.0,
        "n_pairs_used": len(pair_rates),
    }
    return {
        "velocity_mm_per_year": median_rate,
        "n_valid_pairs_per_pixel": n_valid,
        "stats": stats,
    }


def coherence_decay_model(stack: Dict) -> Dict:
    """
    相干性随时间基线衰减模型 — 用全图均值拟合 exp(-baseline/tau)。

    返回 tau(去相干特征时间,天):tau 越大表示场景越稳定。
    """
    if stack["n_pairs"] == 0:
        return {"tau_days": None, "samples": 0}

    xs = []  # 时间基线
    ys = []  # 平均相干性
    for meta, coh in zip(stack["metas"], stack["coherences"]):
        if coh is None:
            continue
        bl = meta.get("temporal_baseline_days")
        if not bl or bl <= 0:
            continue
        m = float(np.nanmean(coh))
        if np.isfinite(m) and m > 0:
            xs.append(bl)
            ys.append(m)

    if len(xs) < 3:
        return {"tau_days": None, "samples": len(xs), "note": "样本不足以拟合"}

    # 拟合 y = exp(-x/tau) → ln(y) = -x/tau
    xs = np.asarray(xs, dtype=np.float32)
    ys = np.asarray(ys, dtype=np.float32)
    ys = np.clip(ys, 1e-6, 1.0)
    lny = np.log(ys)
    # 1/tau ≈ -mean(lny)/mean(x) 近似(更严格应该用 least squares)
    slope, _ = np.polyfit(xs, lny, 1)
    tau = -1.0 / slope if slope < 0 else None
    return {
        "tau_days": float(tau) if tau is not None else None,
        "samples": len(xs),
        "raw_baselines_days": xs.tolist(),
        "raw_mean_coherence": ys.tolist(),
    }
