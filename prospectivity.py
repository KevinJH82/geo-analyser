"""
prospectivity.py — 多证据成矿预测融合(丛林模式 capstone)

在密林/丛林覆盖区,直接蚀变(ratio/pca/band_depth)往往近空白;本模块把多个归一化
证据层加权融合为"成矿有利度"(mineral prospectivity),即便蚀变得分为空也能输出靶区。

证据层(任意尺度,内部统一 percentile 截断归一到 [0,1]):
  E1 alteration   蚀变综合得分 score_map(常规光学)
  E2 veg_stress   地植物学胁迫(红边蓝移/叶绿素↓)—— 丛林区主力证据
  E3 structure    构造邻近度 / 交汇点密度(geo-stru,断裂控矿)
  E4 deformation  InSAR 形变证据 / 沿构造低相干(geo-insar)
  E5 terrain      地形/汇水(可选)

设计原则:
  - 任何证据层缺失(None / 全 NaN / 常数)自动忽略,其权重在"实际存在的层"上再归一(和=1),
    因此丛林区即使只有 veg_stress + structure 也能给出合理有利度。
  - 与现有 insar_analysis.fusion_deformation_mineral 的 _normalize 范式保持一致。
  - 纯函数 + dataclass 返回,便于单测与落盘(沿用 AlterationResult 风格)。
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


# 各证据层默认权重(可被调用方/UI 覆盖)。仅对"实际提供"的层生效并再归一。
DEFAULT_WEIGHTS: Dict[str, float] = {
    "alteration":  0.30,
    "veg_stress":  0.25,
    "structure":   0.20,
    "deformation": 0.15,
    "terrain":     0.10,
}


def normalize_layer(arr: np.ndarray, p_lo: float = 2.0, p_hi: float = 98.0) -> np.ndarray:
    """
    percentile 截断归一到 [0,1]。全 NaN / 常数层返回全 0(等价于"无信息")。
    NaN 像元归一后置 0,以便加权求和;ROI 限定由调用方在融合后统一施加。
    """
    a = np.asarray(arr, dtype=np.float32)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros(a.shape, dtype=np.float32)
    lo = float(np.percentile(finite, p_lo))
    hi = float(np.percentile(finite, p_hi))
    if hi - lo < 1e-9:
        return np.zeros(a.shape, dtype=np.float32)
    out = np.clip((a - lo) / (hi - lo), 0.0, 1.0)
    return np.where(np.isfinite(out), out, 0.0).astype(np.float32)


def _is_informative(arr: Optional[np.ndarray]) -> bool:
    """层是否含信息:非 None、形状有效、存在有限值且非全等。"""
    if arr is None:
        return False
    a = np.asarray(arr, dtype=np.float32)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return False
    return float(finite.max() - finite.min()) > 1e-9


@dataclass
class ProspectivityResult:
    """成矿有利度融合结果(与 InsarAnalysisResult / AlterationResult 风格一致)。"""
    score: np.ndarray                                   # (H,W) [0,1],ROI 外 NaN
    mask: np.ndarray                                    # 高分靶区(top 分位)布尔图
    threshold: float                                    # 高分阈值
    contributions: Dict[str, float] = field(default_factory=dict)  # 各层实际生效权重
    used_layers: List[str] = field(default_factory=list)
    stats: Dict = field(default_factory=dict)
    colormap: str = "inferno"
    unit: str = "成矿有利度 (0-1)"


def fuse_evidence(
    layers: Dict[str, Optional[np.ndarray]],
    weights: Optional[Dict[str, float]] = None,
    roi_mask: Optional[np.ndarray] = None,
    high_percentile: float = 90.0,
    method: str = "weighted_sum",
    gamma: float = 0.9,
) -> ProspectivityResult:
    """
    把多个证据层融合为成矿有利度。

    Parameters
    ----------
    layers : {名称: 2D 数组 | None}
        证据层字典。None / 全空 / 常数层自动忽略。
    weights : {名称: float}, optional
        各层权重;缺省用 DEFAULT_WEIGHTS。仅对实际生效的层再归一(和=1)。
        未在权重表中的层名按 0.1 兜底权重。
    roi_mask : 2D bool, optional
        分析区掩膜;融合后 ROI 外置 NaN。
    high_percentile : float
        高分靶区阈值分位(ROI 内有限像元),默认 top 10%。
    method : {"weighted_sum", "fuzzy_gamma"}
        weighted_sum: 线性加权和(默认,直观稳健)。
        fuzzy_gamma : 模糊 γ 算子 = (fuzzy_sum)^γ · (fuzzy_prod)^(1-γ),
                      对"多证据同时高"更敏感(找矿常用)。

    Returns
    -------
    ProspectivityResult
    """
    base_w = dict(DEFAULT_WEIGHTS)
    if weights:
        base_w.update(weights)

    # 1. 过滤出有信息的层并归一
    present: Dict[str, np.ndarray] = {}
    shape = None
    for name, arr in layers.items():
        if _is_informative(arr):
            present[name] = normalize_layer(arr)
            shape = present[name].shape
    if not present:
        # 全部缺失:返回零有利度(不报错,交由调用方决定)
        h, w = (roi_mask.shape if roi_mask is not None else (1, 1))
        zero = np.zeros((h, w), dtype=np.float32)
        return ProspectivityResult(
            score=np.where(roi_mask, 0.0, np.nan).astype(np.float32) if roi_mask is not None else zero,
            mask=np.zeros((h, w), dtype=bool), threshold=float("nan"),
            contributions={}, used_layers=[],
            stats={"warning": "无可用证据层"},
        )

    # 2. 权重再归一(仅在实际存在的层上)
    eff_w = {name: float(base_w.get(name, 0.1)) for name in present}
    wsum = sum(eff_w.values()) or 1.0
    eff_w = {k: v / wsum for k, v in eff_w.items()}

    # 3. 融合
    if method == "fuzzy_gamma":
        fuzzy_sum = np.zeros(shape, dtype=np.float32)   # 1 - Π(1 - xi)
        prod = np.ones(shape, dtype=np.float32)         # Π xi
        acc = np.ones(shape, dtype=np.float32)
        for name, n in present.items():
            acc *= (1.0 - n)
            prod *= n
        fuzzy_sum = 1.0 - acc
        score = (np.power(np.clip(fuzzy_sum, 0, 1), gamma) *
                 np.power(np.clip(prod, 0, 1), 1.0 - gamma)).astype(np.float32)
    else:  # weighted_sum
        score = np.zeros(shape, dtype=np.float32)
        for name, n in present.items():
            score += eff_w[name] * n
        score = np.clip(score, 0.0, 1.0).astype(np.float32)

    # 4. ROI 限定 + 高分靶区
    if roi_mask is not None:
        score = np.where(roi_mask, score, np.nan).astype(np.float32)
        sample = score[roi_mask & np.isfinite(score)]
    else:
        sample = score[np.isfinite(score)]

    thr = float(np.percentile(sample, high_percentile)) if sample.size else float("nan")
    mask = np.isfinite(score) & (score >= thr) if np.isfinite(thr) else np.zeros(score.shape, dtype=bool)
    if roi_mask is not None:
        mask &= roi_mask

    stats = {
        "n_layers": len(present),
        "high_threshold": thr,
        "high_pixel_ratio": float(mask.sum()) / max(int(sample.size), 1),
        "method": method,
        "score_mean": float(np.nanmean(score)) if sample.size else 0.0,
    }
    return ProspectivityResult(
        score=score, mask=mask, threshold=thr,
        contributions=eff_w, used_layers=list(present.keys()), stats=stats,
    )


# ─────────────────────────────────────────────
# 演示 / 单元测试
# ─────────────────────────────────────────────

def _demo():
    np.random.seed(7)
    H, W = 80, 80
    roi = np.ones((H, W), dtype=bool)

    # 蚀变近空(丛林):全噪声,无结构
    alteration = np.random.uniform(0, 0.05, (H, W)).astype(np.float32)
    # 地植物学胁迫:在 (30:45,30:45) 有强胁迫
    veg = np.random.uniform(0, 0.1, (H, W)).astype(np.float32)
    veg[30:45, 30:45] = 0.9
    # 构造邻近度:在 (35:50,28:42) 有断裂带
    struct = np.zeros((H, W), dtype=np.float32)
    struct[35:50, 28:42] = 0.8
    # InSAR 形变:缺失
    deformation = None

    res = fuse_evidence(
        {"alteration": alteration, "veg_stress": veg,
         "structure": struct, "deformation": deformation},
        roi_mask=roi, high_percentile=92.0, method="weighted_sum",
    )
    print("used_layers:", res.used_layers)
    print("contributions:", {k: round(v, 3) for k, v in res.contributions.items()})
    print("high_threshold: %.3f  high_ratio: %.3f" % (res.threshold, res.stats["high_pixel_ratio"]))

    # veg + structure 交叠区 (35:45,30:42) 应是高有利度靶区
    overlap = res.mask[35:45, 30:42].mean()
    corner = res.mask[0:15, 0:15].mean()
    print("交叠靶区命中率: %.3f   背景角落命中率: %.3f" % (overlap, corner))
    assert "deformation" not in res.used_layers, "缺失层应被忽略"
    assert abs(sum(res.contributions.values()) - 1.0) < 1e-5, "权重应再归一为 1"
    assert overlap > 0.5, "多证据交叠区应判为靶区"
    assert corner < 0.05, "背景不应判为靶区"

    # fuzzy_gamma 对"多证据同时高"更敏感
    res2 = fuse_evidence(
        {"veg_stress": veg, "structure": struct}, roi_mask=roi, method="fuzzy_gamma")
    print("fuzzy_gamma 交叠靶区命中率: %.3f" % res2.mask[35:45, 30:42].mean())

    # 全空层:不报错,零有利度
    res3 = fuse_evidence({"alteration": None}, roi_mask=roi)
    assert res3.used_layers == [] and res3.mask.sum() == 0
    print("PASS: prospectivity.fuse_evidence 工作正常")


if __name__ == "__main__":
    _demo()
