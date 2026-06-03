"""
矿床类型自动识别

输入: ROI 的 GeoJSON 几何
输出: top-k 候选矿床类型,每个含 confidence + evidence + source

两段式判断:
  1. USGS MRDS WFS API — 查 ROI 周边已知矿床点位
  2. Claude API (claude-opus-4-7) — 基于地质背景文字推理(JSON 结构化输出)

合并策略: 数据库权重 0.6 + LLM 权重 0.4,按 deposit_type 名归并后排序。
失败/无 key/无网络时优雅降级,前端再走手动两级下拉。
"""

from __future__ import annotations

import json
import logging
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from alteration_db import all_deposit_type_names

logger = logging.getLogger(__name__)

USGS_MRDS_ENDPOINT = os.environ.get(
    "USGS_MRDS_ENDPOINT",
    "https://mrdata.usgs.gov/services/wfs/mrds",
)
USGS_BUFFER_DEG = 0.5         # ROI bbox 向外扩展(度,约 50km)
USGS_TIMEOUT = 10             # 秒
USGS_MAX_FEATURES = 200

DEEPSEEK_MODEL = os.environ.get("DEEPSEEK_MODEL", "deepseek-chat")
DEEPSEEK_BASE_URL = os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com")
LLM_TIMEOUT = 30
LLM_MAX_TOKENS = 1500

# 数据库 / LLM 合并权重
W_DB = 0.6
W_LLM = 0.4


# ─────────────────────────────────────────────
# 几何工具
# ─────────────────────────────────────────────

def _roi_bbox(roi_geojson: Dict[str, Any]) -> Optional[Tuple[float, float, float, float]]:
    """从 GeoJSON 几何提取 bbox (minLon, minLat, maxLon, maxLat)。"""
    geom = roi_geojson.get("geometry") if "geometry" in roi_geojson else roi_geojson
    if not geom:
        return None
    gtype = geom.get("type")
    coords = geom.get("coordinates")
    if not coords:
        return None

    pts: List[Tuple[float, float]] = []
    if gtype == "Polygon":
        for ring in coords:
            pts.extend((float(x), float(y)) for x, y in ring)
    elif gtype == "MultiPolygon":
        for poly in coords:
            for ring in poly:
                pts.extend((float(x), float(y)) for x, y in ring)
    elif gtype == "Point":
        pts.append((float(coords[0]), float(coords[1])))
    else:
        return None

    if not pts:
        return None
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    return min(lons), min(lats), max(lons), max(lats)


# ─────────────────────────────────────────────
# 第一段: USGS MRDS 查询
# ─────────────────────────────────────────────

# USGS MRDS dep_type / commodity 到本地 JSON 中 deposit_type 名称的关键词映射
# 用模糊匹配 (子串 / 关键词组),不要求严格相等。
_USGS_KEYWORDS = [
    # (USGS 字段中的关键词列表, 本地 JSON 中的 deposit_type_name)
    (["porphyry", "porphyry cu", "斑岩"],                    "斑岩型铜矿"),
    (["porphyry mo", "斑岩钼"],                              "斑岩型钼矿（Climax型）"),
    (["porphyry cu-mo", "porphyry copper-molybdenum"],       "斑岩型铜钼矿"),
    (["porphyry cu-au", "porphyry copper-gold"],             "斑岩型铜金矿"),
    (["skarn cu", "skarn copper", "矽卡岩铜"],               "矽卡岩型铜矿"),
    (["skarn fe", "skarn iron", "矽卡岩铁"],                 "矽卡岩型铁矿"),
    (["skarn w", "skarn tungsten", "矽卡岩钨"],              "矽卡岩型钨矿"),
    (["skarn pb-zn", "skarn lead-zinc"],                     "矽卡岩型铅锌矿"),
    (["vms", "volcanogenic massive sulfide"],                "VMS型铜锌矿"),
    (["iocg", "iron oxide copper-gold"],                     "IOCG型铁氧化物铜金矿"),
    (["mvt", "mississippi valley"],                          "MVT型铅锌矿"),
    (["sedex", "sedimentary exhalative"],                    "SEDEX型铅锌矿"),
    (["epithermal", "low-sulfidation", "浅成低温"],          "浅成低温热液型金矿（低硫化型）"),
    (["high-sulfidation", "high sulfidation"],               "浅成低温热液型金矿（高硫化型）"),
    (["orogenic gold", "造山金"],                            "造山型金矿"),
    (["carlin"],                                             "卡林型金矿"),
    (["pge", "platinum-group", "layered intrusion"],         "层状镁铁质侵入体PGE矿"),
    (["bif", "banded iron formation", "条带状铁"],           "BIF型铁矿（条带状铁建造）"),
    (["greisen w-sn", "greisen tungsten-tin"],               "云英岩型钨锡矿"),
    (["greisen sn", "greisen tin"],                          "云英岩型锡矿"),
    (["magmatic ni-cu", "magmatic nickel-copper"],           "岩浆硫化物型镍铜矿"),
    (["laterite ni", "lateritic nickel"],                    "红土型镍矿"),
    (["sediment-hosted cu", "sediment hosted copper"],       "沉积岩容矿铜钴矿"),
    (["laterite co", "lateritic cobalt"],                    "红土型钴矿"),
    (["sedimentary mn", "sedimentary manganese"],            "沉积型锰矿"),
    (["volcanogenic mn", "volcanic-sedimentary mn"],         "火山沉积型锰矿"),
    (["podiform chromite", "ophiolite chromite"],            "豆荚状铬铁矿（蛇绿岩型）"),
    (["carbonatite ree", "carbonatite rare earth"],          "碳酸岩型稀土矿"),
    (["ion adsorption ree"],                                 "离子吸附型稀土矿"),
    (["pegmatite li", "spodumene"],                          "伟晶岩型锂矿（锂辉石型）"),
    (["brine li", "salar"],                                  "盐湖卤水型锂矿"),
    (["clay li", "lithium clay"],                            "沉积型锂矿（粘土型）"),
    (["unconformity uranium"],                               "不整合面型铀矿"),
    (["sandstone uranium"],                                  "砂岩型铀矿"),
    (["albitite uranium", "na-metasomatic"],                 "钠交代型铀矿"),
    (["laterite bauxite"],                                   "红土型铝土矿"),
    (["epithermal sb", "low-temperature antimony"],          "低温热液型锑矿"),
    (["epithermal hg", "low-temperature mercury"],           "低温热液型汞矿"),
    (["magmatic v-ti-fe", "vanadium-titanium-iron"],         "岩浆型钒钛磁铁矿"),
    (["black shale v", "shale vanadium"],                    "沉积型钒矿（黑色页岩型）"),
    (["graphite"],                                           "区域变质型石墨矿"),
    (["kimberlite diamond"],                                 "金伯利岩型金刚石矿"),
    (["evaporite potash"],                                   "蒸发岩型钾盐矿"),
    (["phosphate"],                                          "沉积型磷块岩矿"),
    (["fluorite vein"],                                      "热液脉型萤石矿"),
    (["barite"],                                             "热液/沉积型重晶石矿"),
]


def _match_usgs_record(commodity: str, dep_type: str) -> Optional[str]:
    """把 USGS 的 commodity + dep_type 字段映射到本地 deposit_type 名称。"""
    blob = f"{commodity} {dep_type}".lower()
    for keywords, name in _USGS_KEYWORDS:
        for kw in keywords:
            if kw.lower() in blob:
                return name
    return None


def _query_usgs_mrds(bbox: Tuple[float, float, float, float]) -> List[Dict[str, str]]:
    """
    查询 USGS MRDS WFS 返回 bbox 内的矿床点。
    返回 [{commodity, dep_type, name}, ...] 列表。
    """
    try:
        import requests
    except ImportError:
        logger.warning("requests 未安装,无法查询 USGS")
        return []

    min_lon, min_lat, max_lon, max_lat = bbox
    # 向外扩展 buffer
    min_lon -= USGS_BUFFER_DEG
    max_lon += USGS_BUFFER_DEG
    min_lat -= USGS_BUFFER_DEG
    max_lat += USGS_BUFFER_DEG

    params = {
        "service":      "WFS",
        "version":      "1.1.0",
        "request":      "GetFeature",
        "typeName":     "mrds",
        "outputFormat": "application/json",
        "bbox":         f"{min_lon},{min_lat},{max_lon},{max_lat},EPSG:4326",
        "maxFeatures":  USGS_MAX_FEATURES,
    }
    try:
        r = requests.get(USGS_MRDS_ENDPOINT, params=params, timeout=USGS_TIMEOUT)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        logger.warning(f"USGS MRDS 查询失败: {e}")
        return []

    out = []
    for feat in data.get("features", []):
        props = feat.get("properties", {}) or {}
        commodity = str(props.get("commod1", "") or "")
        dep_type = str(props.get("dep_type", "") or "")
        name = str(props.get("site_name", "") or "")
        if commodity or dep_type:
            out.append({"commodity": commodity, "dep_type": dep_type, "name": name})
    return out


def _infer_from_usgs(roi_geojson: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    返回 {deposit_type_name: {confidence, evidence, source, hit_count}}
    confidence: 命中数越多越高,clipped to [0, 1]
    """
    bbox = _roi_bbox(roi_geojson)
    if not bbox:
        return {}

    records = _query_usgs_mrds(bbox)
    if not records:
        return {}

    # 命中聚合
    type_hits: Dict[str, List[str]] = defaultdict(list)
    for r in records:
        name = _match_usgs_record(r["commodity"], r["dep_type"])
        if name:
            label = r["name"] or f"{r['commodity']} {r['dep_type']}".strip()
            type_hits[name].append(label)

    out = {}
    for name, hits in type_hits.items():
        # confidence: 1 个命中 0.55,3 个 0.80,5+ 个 0.95
        n = len(hits)
        if n == 1:
            conf = 0.55
        elif n == 2:
            conf = 0.70
        elif n <= 4:
            conf = 0.85
        else:
            conf = 0.95
        sample_names = ", ".join(hits[:3])
        evidence = f"USGS MRDS 在 ROI 周边({USGS_BUFFER_DEG}°缓冲)发现 {n} 个相关矿点(示例: {sample_names})"
        out[name] = {
            "confidence": conf, "evidence": evidence,
            "source": "USGS MRDS", "hit_count": n,
        }
    return out


# ─────────────────────────────────────────────
# 第二段: Claude API 推理
# ─────────────────────────────────────────────

_LLM_SYSTEM_TEMPLATE = """你是资深地质矿产专家,擅长根据地理位置和地质背景判断该区域可能发育的矿床类型。

任务: 给定一个 ROI 区域(经纬度边界),从下面给定的矿床类型枚举中,推荐该区域最有可能发育的 top-{top_k} 个矿床类型。

【优先级 1 — 已知成矿带/矿集区直接命中】
如果 ROI 落在以下任一已知成矿带/著名矿集区内,**该带的代表矿床类型必须排在第一**,confidence 不低于 0.85,evidence 中明确点出矿带名称:

中国典型成矿带:
- 招远-莱州金矿带 / 玲珑金矿带 / 焦家金矿带 (山东半岛胶东) → 造山型金矿
- 长江中下游铜铁多金属矿带 (湖北鄂东南、安徽铜陵、江苏宁芜) → 矽卡岩型铜矿/矽卡岩型铁矿/斑岩型铜钼矿
- 冈底斯斑岩铜矿带 (西藏中南部) → 斑岩型铜矿
- 玉龙斑岩铜矿带 (青藏三江) → 斑岩型铜矿
- 钦杭成矿带 (湖南南岭、江西) → 云英岩型钨锡矿/矽卡岩型钨矿/伟晶岩型锂矿
- 三江多金属成矿带 (川滇藏) → 浅成低温热液型金矿/斑岩型铜矿
- 东天山-北山多金属带 (新疆) → 斑岩型铜矿/VMS型铜锌矿
- 鞍山-本溪铁矿带 (辽宁) → BIF型铁矿
- 攀西V-Ti磁铁矿带 (四川攀枝花) → 岩浆型钒钛磁铁矿
- 大兴安岭多金属带 (内蒙古东北) → 斑岩型钼矿/浅成低温热液型银矿
- 滇黔桂金三角 (贵州西南、云南东南) → 卡林型金矿
- 鄂尔多斯盆地砂岩铀矿带 (内蒙古、陕西、新疆) → 砂岩型铀矿
- 柴达木盐湖矿集区 (青海) → 盐湖卤水型锂矿/蒸发岩型钾盐矿

国际典型成矿带(同样规则):
- 智利-秘鲁安第斯斑岩铜矿带 → 斑岩型铜矿
- 美国西部造山带 → 卡林型金矿、斑岩型铜钼矿
- 加拿大地盾 → BIF型铁矿、VMS型铜锌矿
- 西伯利亚地台 → 岩浆硫化物型镍铜矿、金伯利岩型金刚石矿

【优先级 2 — 大地构造 + 地质背景】
ROI 不在已知矿带时,基于:
- 大地构造背景(克拉通/造山带/岛弧/裂谷/被动陆缘)
- 区域地质(岩浆活动、变质程度、沉积环境)
推断最可能的矿床类型,confidence 0.4~0.7。

【约束】
矿床类型必须严格从以下列表中选择(不要发明,不要列表外名称):
{deposit_type_enum}

对每个推荐:
- confidence: 0~1,反映该区域发育该矿床类型的可能性
- evidence: 1~2 句话,引用具体地质事实和矿带名,如"位于胶东招远-莱州金矿带,中生代花岗岩中受 NE 向断裂控制"
"""

# OpenAI function-calling 格式(DeepSeek 兼容)
_LLM_TOOL_SCHEMA = {
    "type": "function",
    "function": {
        "name": "report_deposit_types",
        "description": "返回 top-k 候选矿床类型,按 confidence 降序。",
        "parameters": {
            "type": "object",
            "properties": {
                "candidates": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "deposit_type": {"type": "string"},
                            "confidence":   {"type": "number", "minimum": 0, "maximum": 1},
                            "evidence":     {"type": "string"},
                        },
                        "required": ["deposit_type", "confidence", "evidence"],
                    },
                },
            },
            "required": ["candidates"],
        },
    },
}


def _infer_from_llm(roi_geojson: Dict[str, Any], top_k: int = 3) -> Dict[str, Dict[str, Any]]:
    """
    调 DeepSeek (deepseek-chat) 推理,OpenAI 兼容 API。
    返回 {deposit_type_name: {confidence, evidence, source}}。
    无 API key / SDK / 网络失败时返回 {}。
    """
    api_key = os.environ.get("DEEPSEEK_API_KEY", "").strip()
    if not api_key:
        logger.info("DEEPSEEK_API_KEY 未设置,跳过 LLM 推理")
        return {}

    try:
        from openai import OpenAI
    except ImportError:
        logger.warning("openai SDK 未安装,跳过 LLM 推理")
        return {}

    bbox = _roi_bbox(roi_geojson)
    if not bbox:
        return {}
    min_lon, min_lat, max_lon, max_lat = bbox
    center_lon = (min_lon + max_lon) / 2
    center_lat = (min_lat + max_lat) / 2

    enum_str = "\n".join(f"- {n}" for n in all_deposit_type_names())
    system_text = _LLM_SYSTEM_TEMPLATE.format(top_k=top_k, deposit_type_enum=enum_str)

    user_text = (
        f"ROI 边界: 经度 [{min_lon:.4f}, {max_lon:.4f}],"
        f"纬度 [{min_lat:.4f}, {max_lat:.4f}]\n"
        f"ROI 中心: ({center_lon:.4f}°E, {center_lat:.4f}°N)\n\n"
        f"请调用 report_deposit_types 工具,给出 top-{top_k} 候选。"
    )

    try:
        client = OpenAI(api_key=api_key, base_url=DEEPSEEK_BASE_URL, timeout=LLM_TIMEOUT)
        resp = client.chat.completions.create(
            model=DEEPSEEK_MODEL,
            max_tokens=LLM_MAX_TOKENS,
            messages=[
                {"role": "system", "content": system_text},
                {"role": "user",   "content": user_text},
            ],
            tools=[_LLM_TOOL_SCHEMA],
            tool_choice={"type": "function", "function": {"name": "report_deposit_types"}},
        )
    except Exception as e:
        logger.warning(f"DeepSeek API 调用失败: {e}")
        return {}

    # 提取 tool_call 输出
    candidates = []
    try:
        msg = resp.choices[0].message
        for tc in (msg.tool_calls or []):
            if tc.function and tc.function.name == "report_deposit_types":
                args = json.loads(tc.function.arguments or "{}")
                candidates = args.get("candidates", [])
                break
    except Exception as e:
        logger.warning(f"DeepSeek 响应解析失败: {e}")
        return {}

    valid_names = set(all_deposit_type_names())
    out = {}
    for c in candidates[:top_k]:
        name = c.get("deposit_type", "").strip()
        if name not in valid_names:
            continue
        conf = float(c.get("confidence", 0.5))
        conf = max(0.0, min(1.0, conf))
        evidence = str(c.get("evidence", "")).strip()
        out[name] = {
            "confidence": conf, "evidence": evidence,
            "source": f"DeepSeek ({DEEPSEEK_MODEL})",
        }
    return out


# ─────────────────────────────────────────────
# 合并 + 主入口
# ─────────────────────────────────────────────

def infer_deposit_types(roi_geojson: Dict[str, Any], top_k: int = 3) -> List[Dict[str, Any]]:
    """
    主入口: 返回 top-k 候选矿床类型列表(按 confidence 降序)。
    每项: {deposit_type, confidence, evidence, source}
    """
    db_results = _infer_from_usgs(roi_geojson)
    llm_results = _infer_from_llm(roi_geojson, top_k=max(top_k, 3))

    # 动态权重:
    # - 两边都有 → 0.6 × DB + 0.4 × LLM(数据库为主)
    # - 只有 LLM(USGS 无返回/网络失败) → LLM 权重提到 1.0,直接用原始置信度
    # - 只有 DB → 用 DB 原始置信度(不再×0.6,避免无谓衰减)
    merged: Dict[str, Dict[str, Any]] = {}
    all_names = set(db_results) | set(llm_results)
    for name in all_names:
        db = db_results.get(name)
        llm = llm_results.get(name)
        if db and llm:
            conf = W_DB * db["confidence"] + W_LLM * llm["confidence"]
            evidence = f"{db['evidence']};同时 AI 推理 {llm['evidence']}"
            source = f"{db['source']} + {llm['source']}"
        elif db:
            conf = db["confidence"]
            evidence = db["evidence"]
            source = db["source"]
        else:  # LLM only
            conf = llm["confidence"]
            evidence = llm["evidence"]
            source = llm["source"]
        merged[name] = {
            "deposit_type": name,
            "confidence":   round(conf, 3),
            "evidence":     evidence,
            "source":       source,
        }

    sorted_list = sorted(merged.values(), key=lambda x: -x["confidence"])
    return sorted_list[:top_k]


if __name__ == "__main__":
    # 演示: 智利斑岩铜矿带某处 (Chuquicamata 附近)
    demo_roi = {
        "type": "Polygon",
        "coordinates": [[
            [-68.95, -22.30], [-68.85, -22.30],
            [-68.85, -22.40], [-68.95, -22.40],
            [-68.95, -22.30],
        ]],
    }
    print("=== Chuquicamata (智利斑岩铜矿带) ===")
    for c in infer_deposit_types(demo_roi, top_k=3):
        print(f"  [{c['confidence']:.2f}] {c['deposit_type']}  ({c['source']})")
        print(f"         {c['evidence']}")
