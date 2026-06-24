"""
structural_deposit.py — 从 geo-stru 构造解译产出获取矿床类型(最高优先级来源)

geo-stru 已基于断裂走向/密度/地形/形变做了纯构造的矿床类型推理,落盘在
<AOI>/structural/<run_id>/metadata.json 的 deposit_inference 字段,经
commons/structural_broker 按 bbox 相交发现(与 structural_weighting.py 同一订阅链路)。

本模块只读消费这份产出,并把 geo-stru 的「构造控矿类型」翻译成 geo-analyser 蚀变库的
type_name —— 二者命名体系不一致(geo-stru 10 个粗粒度构造控矿类型,蚀变库 ~60 个细分类型),
故用一张显式映射表(GEO_STRU_TO_DB_TYPES,键对齐 geo-stru/core/deposit_inference.py 的
DEPOSIT_RULES)。油气/煤层气类型映射到蚀变库的「微渗漏蚀变」类型(烃微渗漏在近地表的
红层褪色/次生碳酸盐化/粘土化等遥感异常);仍无对应的类型显式留空 → 丢弃。

输出形态刻意与 deposit_type_inference.infer_deposit_types_detailed() 的候选一致
({deposit_type, confidence, evidence, source}),供前端 renderDepositCandidates 直接复用。

无 geo-stru 产物 / 全部映射不到时优雅降级(返回空候选 + reason),前端再回退 LLM / 手动。
"""

from __future__ import annotations

import logging
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

from alteration_db import all_deposit_type_names, deposit_type_index
from delivery_project import bbox_from_geojson

logger = logging.getLogger(__name__)

SOURCE_LABEL = "geo-stru构造推理"


# ─────────────────────────────────────────────
# geo-stru 构造控矿类型 → geo-analyser 蚀变库 type_name(一对多,顺序=优先)
# 键须与 geo-stru/core/deposit_inference.py 的 DEPOSIT_RULES 保持同步。
# 值在模块加载时对照 all_deposit_type_names() 校验,库中不存在的名字会被剔除并告警。
# ─────────────────────────────────────────────
GEO_STRU_TO_DB_TYPES: Dict[str, List[str]] = {
    "蚀变岩型金矿(破碎带)":      ["破碎带蚀变岩型金矿（焦家式）", "造山型金矿"],
    "斑岩型铜钼矿":              ["斑岩型铜钼矿", "斑岩型铜矿", "斑岩型钼矿（Climax型）", "斑岩型铜金矿"],
    "矽卡岩型铁矿/铜矿":          ["矽卡岩型铁矿", "矽卡岩型铜矿"],
    "VMS型多金属矿":            ["VMS型铜锌矿"],
    "石英脉型金矿":              ["热液脉型金矿（石英脉型，玲珑式）"],
    "沉积型矿产(煤/铝土/盐类)":   ["红土型铝土矿", "蒸发岩型钾盐矿"],
    "BIF型铁矿":               ["BIF型铁矿（条带状铁建造）"],
    # 油气/煤层气:映射到蚀变库的「微渗漏蚀变」类型 —— 烃类沿微裂隙向上渗漏在近地表
    # 引起红层褪色/次生碳酸盐化/粘土化/植被胁迫等遥感可测异常(见 alteration_db 微渗漏靶点)。
    "常规油气藏(微渗漏)":  ["常规油藏(微渗漏蚀变模式)", "常规气藏(微渗漏蚀变模式)"],
    "致密油气/页岩气":     ["致密油藏(微渗漏蚀变)"],
    "煤层气/煤炭":         ["煤层气(微渗漏蚀变)"],
}


def _validated_mapping() -> Dict[str, List[str]]:
    """加载期校验:剔除蚀变库中不存在的 type_name(防全角括号/改名踩坑)。"""
    valid = set(all_deposit_type_names())
    out: Dict[str, List[str]] = {}
    for stru_type, db_names in GEO_STRU_TO_DB_TYPES.items():
        kept = [n for n in db_names if n in valid]
        dropped = [n for n in db_names if n not in valid]
        if dropped:
            logger.warning(
                "structural_deposit 映射表中以下 type_name 不在蚀变库,已剔除: %s (geo-stru类型=%s)",
                dropped, stru_type,
            )
        out[stru_type] = kept
    return out


# 模块加载时校验一次(失败不阻断导入,退回原表)
try:
    _MAPPING = _validated_mapping()
except Exception as e:  # noqa: BLE001
    logger.warning("structural_deposit 映射表校验失败,使用未校验表: %s", e)
    _MAPPING = dict(GEO_STRU_TO_DB_TYPES)


# ─────────────────────────────────────────────
# geo-stru 产物发现(复用 commons/structural_broker)
# ─────────────────────────────────────────────

def _find_structural(bbox: Tuple[float, float, float, float]) -> List[Dict[str, Any]]:
    """与 structural_weighting 一致:经 commons.structural_broker 按 bbox 相交发现 geo-stru 产物。"""
    if "/opt/deepexplor-services" not in sys.path:
        sys.path.insert(0, "/opt/deepexplor-services")
    try:
        from commons.structural_broker import find_structural_for_bbox
    except Exception as e:  # noqa: BLE001
        logger.info("commons.structural_broker 不可用,跳过 geo-stru 矿床类型: %s", e)
        return []
    try:
        return find_structural_for_bbox(bbox) or []
    except Exception as e:  # noqa: BLE001
        logger.warning("find_structural_for_bbox 失败: %s", e)
        return []


def _overlap_area(a, b) -> float:
    """两个 [min_lon,min_lat,max_lon,max_lat] 的相交面积(度²),不相交返回 0。"""
    if not a or not b or len(a) < 4 or len(b) < 4:
        return 0.0
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    if dx <= 0 or dy <= 0:
        return 0.0
    return dx * dy


def _pick_entry(matches: List[Dict[str, Any]],
                bbox: Tuple[float, float, float, float]) -> Optional[Dict[str, Any]]:
    """多个相交产物时:优先有 deposit_inference 的,再按与 ROI 重叠面积最大。"""
    usable = [m for m in matches if (m.get("deposit_inference") or {}).get("candidates")]
    if not usable:
        return None
    return max(usable, key=lambda m: _overlap_area(m.get("aoi_bbox"), bbox))


# ─────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────

# 收集首选之外的次选候选时,与首选置信度的最大允许差(避免把 geo-stru 给所有类型的
# 基础分误当成"候选";只保留与首选接近的同档类型)。
PRIMARY_CONF_BAND = 0.12


def _evidence_text(sc: Dict[str, Any], fallback: str = "") -> str:
    ev_list = sc.get("evidence") or []
    ev = "；".join(str(x) for x in ev_list) if isinstance(ev_list, list) else str(ev_list)
    control_model = sc.get("control_model") or ""
    if control_model:
        ev = f"{ev}｜控矿模式:{control_model}" if ev else f"控矿模式:{control_model}"
    return ev or fallback


def structural_deposit_candidates(roi_geojson: Dict[str, Any], top_k: int = 3) -> Dict[str, Any]:
    """
    从 geo-stru 构造推理产出获取矿床类型候选(翻译成蚀变库 type_name)。

    核心原则:**锚定 geo-stru 的首选(置信度最高)结论** —— 它代表 geo-stru 对该区矿床
    类型的权威判断。geo-stru 会给全部 10 个类型打分(油气区里金属类型也有 0.5~0.6 基础分),
    因此绝不能在首选被丢弃后顺位提升一个低分金属类型(会把"石油"悄悄改判成"铜锌")。

      - 首选可映射到蚀变库 → status="ok",仅收集与首选同档(置信度接近)的可映射候选。
      - 首选不可映射(油气/煤层气/煤炭等非蚀变靶区)→ status="non_alteration",
        不改判、不回退 LLM,明确告知蚀变分析不适用。
      - 无 geo-stru 产物 / 无候选 → status="no_product",前端回退 LLM。

    Returns
    -------
    {
        "candidates": [{deposit_type, confidence, evidence, source}, ...],  # 蚀变库 type_name
        "status":     "ok" | "no_product" | "non_alteration",
        "reason":     str,            # 降级/不适用原因(供前端 hint)
        "source":     "geo-stru",
        "aoi_name":   str|None,
        "run_id":     str|None,
        "primary_model": str|None,    # geo-stru 原始首选类型(展示用)
        "primary_confidence": float|None,
    }
    """
    base = {"candidates": [], "status": "no_product", "reason": "", "source": "geo-stru",
            "aoi_name": None, "run_id": None, "primary_model": None, "primary_confidence": None}

    bbox = bbox_from_geojson(roi_geojson) if roi_geojson else None
    if not bbox:
        base["reason"] = "ROI 几何无效,无法匹配 geo-stru 构造产物"
        return base

    entry = _pick_entry(_find_structural(bbox), bbox)
    if entry is None:
        base["reason"] = "该区无 geo-stru 构造解译产物(含矿床类型推理),回退 LLM 自动识别"
        return base

    di = entry.get("deposit_inference") or {}
    summary = di.get("structural_control_summary") or ""
    stru_cands = di.get("candidates") or []
    base["aoi_name"] = entry.get("aoi_name")
    base["run_id"] = entry.get("run_id")

    if not stru_cands:
        base["reason"] = "geo-stru 构造产物无矿床类型候选,回退 LLM 自动识别"
        return base

    idx = deposit_type_index()

    # 锚定首选(候选已按 confidence 降序,stru_cands[0] 即最高分)
    primary = stru_cands[0]
    primary_type = (primary.get("deposit_type") or "").strip()
    primary_conf = float(primary.get("confidence", 0.0))
    base["primary_model"] = primary_type
    base["primary_confidence"] = round(primary_conf, 3)

    primary_db = [n for n in _MAPPING.get(primary_type, []) if n in idx]
    if not primary_db:
        # 首选为油气/煤层气/煤炭等非蚀变靶区:不改判金属矿,也不回退 LLM
        base["status"] = "non_alteration"
        base["reason"] = (
            f"geo-stru 判定该区首选矿床类型为「{primary_type}」(置信度 {primary_conf:.2f}),"
            f"属油气/非蚀变类靶区,蚀变分析不适用;如确需分析请手动选择矿种与矿床类型。"
        )
        return base

    # 首选为蚀变类:只收集与首选同档(置信度差 ≤ PRIMARY_CONF_BAND)的可映射候选,
    # 去重保留最高 confidence。首选自身的映射因分最高自然排在最前。
    best: Dict[str, Dict[str, Any]] = {}
    for sc in stru_cands:
        conf = float(sc.get("confidence", 0.0))
        if primary_conf - conf > PRIMARY_CONF_BAND:
            continue
        stru_type = (sc.get("deposit_type") or "").strip()
        db_names = [n for n in _MAPPING.get(stru_type, []) if n in idx]
        if not db_names:
            continue
        evidence = _evidence_text(sc, fallback=summary)
        for db_name in db_names:
            prev = best.get(db_name)
            if prev is None or conf > prev["confidence"]:
                best[db_name] = {
                    "deposit_type": db_name,
                    "confidence": round(conf, 3),
                    "evidence": evidence,
                    "source": SOURCE_LABEL,
                }

    candidates = sorted(best.values(), key=lambda c: -c["confidence"])[:top_k]
    base["candidates"] = candidates
    base["status"] = "ok"
    return base


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("=== 映射表校验 ===")
    for k, v in _MAPPING.items():
        print(f"  {k} → {v}")
