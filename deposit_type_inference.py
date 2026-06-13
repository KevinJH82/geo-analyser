"""
矿床类型自动识别

输入: ROI 的 GeoJSON 几何
输出: top-k 候选矿床类型,每个含 confidence + evidence + source

机制:
  1. LLM 推理 (DeepSeek, OpenAI 兼容) — 主路径,基于地理坐标 + 成矿带知识
     从本地库的 deposit_type 枚举中给出 top-k 候选(结构化工具调用输出)。
  2. USGS MRDS WFS — 辅助佐证。该 WFS 图层只提供周边已知矿点的"矿种代码 + 名称"
     (无矿床类型字段),故仅用作:① 喂给 LLM 提示词;② LLM 候选矿种命中周边矿种时
     提升置信度。它本身给不出 deposit_type。

失败/无 key/无网络时优雅降级,返回空候选 + reason,前端再走手动两级下拉。
"""

from __future__ import annotations

import json
import logging
import os
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

from alteration_db import all_deposit_type_names, get_deposit_type_meta

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

# USGS 周边矿种命中时给 LLM 候选的置信度加成(封顶 1.0)
USGS_CORROBORATION_BONUS = 0.1

# USGS MRDS code_list 矿种代码(元素符号)→ 本地库中文矿种名
# 用于把"周边已知矿点"对齐到 deposit_type 的归属矿种,做佐证加权。
_CODE_TO_COMMODITY = {
    "AU": "金", "AG": "银", "CU": "铜", "PB": "铅锌", "ZN": "铅锌",
    "MO": "钼", "W": "钨", "SN": "锡", "FE": "铁", "MN": "锰",
    "TI": "钛", "CR": "铬", "V": "钒", "CO": "钴", "NI": "镍",
    "HG": "汞", "SB": "锑", "U": "铀", "LI": "锂", "TA": "钽",
    "REE": "稀土", "CE": "稀土", "LA": "稀土", "Y": "稀土", "ND": "稀土",
    "P": "磷矿", "C": "石墨", "F": "萤石", "BA": "重晶石", "AL": "铝土矿",
    "PT": "铂族元素", "PD": "铂族元素", "K": "钾盐",
}


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
# 辅助佐证: USGS MRDS 周边矿种
# ─────────────────────────────────────────────
# 注意: USGS MRDS WFS 的 mrds 图层只暴露 site_name / code_list 等字段,
# 没有 dep_type(矿床类型)字段,故只能给出"周边已知矿点的矿种 + 名称",
# 给不出矿床类型。本段仅作 LLM 推理的佐证,不产 deposit_type。

_MS_NS = "{http://mapserver.gis.umn.edu/mapserver}"


def _query_usgs_mrds(bbox: Tuple[float, float, float, float]) -> List[Dict[str, Any]]:
    """
    查询 USGS MRDS WFS 返回 bbox 内的矿床点。
    返回 [{commodity_codes: [..], name}, ...] 列表。

    要点(均经在线服务实测):
      - WFS GetFeature 仅支持 GML 输出(text/xml; subtype=gml/3.1.1),不支持 JSON。
      - version=1.1.0 + EPSG:4326 时 bbox 轴序为 lat,lon(不是 lon,lat)。
      - mrds 图层无 dep_type/commod1 字段,矿种在 code_list(空格分隔的元素代码)。
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
        "outputFormat": "text/xml; subtype=gml/3.1.1",
        # WFS 1.1.0 + EPSG:4326 轴序为 lat,lon
        "bbox":         f"{min_lat},{min_lon},{max_lat},{max_lon},EPSG:4326",
        "maxFeatures":  USGS_MAX_FEATURES,
    }
    try:
        r = requests.get(USGS_MRDS_ENDPOINT, params=params, timeout=USGS_TIMEOUT)
        r.raise_for_status()
        root = ET.fromstring(r.content)
    except Exception as e:
        logger.warning(f"USGS MRDS 查询失败: {e}")
        return []

    out: List[Dict[str, Any]] = []
    for feat in root.iter(f"{_MS_NS}mrds"):
        name_el = feat.find(f"{_MS_NS}site_name")
        code_el = feat.find(f"{_MS_NS}code_list")
        name = (name_el.text or "").strip() if name_el is not None else ""
        codes_raw = (code_el.text or "").strip() if code_el is not None else ""
        codes = [c.strip().upper() for c in codes_raw.split() if c.strip()]
        if codes or name:
            out.append({"commodity_codes": codes, "name": name})
    return out


def _usgs_nearby_context(roi_geojson: Dict[str, Any]) -> Dict[str, Any]:
    """
    返回 ROI 周边已知矿点的矿种佐证:
        {commodities: set[中文矿种名], sample_sites: [矿点名, ...]}
    用于喂给 LLM 提示词 + 对候选做佐证加权。失败/无网络时返回空集合。
    """
    empty = {"commodities": set(), "sample_sites": []}
    bbox = _roi_bbox(roi_geojson)
    if not bbox:
        return empty

    records = _query_usgs_mrds(bbox)
    if not records:
        return empty

    commodities: set = set()
    sample_sites: List[str] = []
    for r in records:
        for code in r["commodity_codes"]:
            cn = _CODE_TO_COMMODITY.get(code)
            if cn:
                commodities.add(cn)
        if r["name"] and r["name"] not in sample_sites and len(sample_sites) < 8:
            sample_sites.append(r["name"])
    return {"commodities": commodities, "sample_sites": sample_sites}


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


def _infer_from_llm(
    roi_geojson: Dict[str, Any],
    top_k: int = 3,
    usgs_context: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, Dict[str, Any]], str]:
    """
    调 DeepSeek (deepseek-chat) 推理,OpenAI 兼容 API。
    返回 (results, reason):
      results = {deposit_type_name: {confidence, evidence, source}}
      reason  = 空串表示成功;否则为降级原因(无 key / 无 SDK / 调用失败等)。
    """
    api_key = os.environ.get("DEEPSEEK_API_KEY", "").strip()
    if not api_key:
        logger.info("DEEPSEEK_API_KEY 未设置,跳过 LLM 推理")
        return {}, "未配置 DEEPSEEK_API_KEY,无法自动识别矿床类型,请手动选择"

    try:
        from openai import OpenAI
    except ImportError:
        logger.warning("openai SDK 未安装,跳过 LLM 推理")
        return {}, "未安装 openai SDK,无法调用 LLM,请手动选择"

    bbox = _roi_bbox(roi_geojson)
    if not bbox:
        return {}, "ROI 几何无效,无法提取坐标"
    min_lon, min_lat, max_lon, max_lat = bbox
    center_lon = (min_lon + max_lon) / 2
    center_lat = (min_lat + max_lat) / 2

    enum_str = "\n".join(f"- {n}" for n in all_deposit_type_names())
    system_text = _LLM_SYSTEM_TEMPLATE.format(top_k=top_k, deposit_type_enum=enum_str)

    user_lines = [
        f"ROI 边界: 经度 [{min_lon:.4f}, {max_lon:.4f}],"
        f"纬度 [{min_lat:.4f}, {max_lat:.4f}]",
        f"ROI 中心: ({center_lon:.4f}°E, {center_lat:.4f}°N)",
    ]
    # USGS 周边矿种佐证(若有)注入提示词
    nearby = (usgs_context or {}).get("commodities") or set()
    sites = (usgs_context or {}).get("sample_sites") or []
    if nearby:
        line = f"参考: USGS MRDS 在 ROI 周边已知矿点涉及矿种: {', '.join(sorted(nearby))}"
        if sites:
            line += f";示例矿点: {', '.join(sites[:5])}"
        line += "。请结合该矿种线索判断,但矿床类型仍须从给定枚举中选择。"
        user_lines.append(line)
    user_lines.append(f"\n请调用 report_deposit_types 工具,给出 top-{top_k} 候选。")
    user_text = "\n".join(user_lines)

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
        return {}, f"LLM 调用失败: {e}"

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
        return {}, f"LLM 响应解析失败: {e}"

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
    return out, ""


# ─────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────

def infer_deposit_types_detailed(roi_geojson: Dict[str, Any], top_k: int = 3) -> Dict[str, Any]:
    """
    主入口(带降级原因): 返回 {candidates, reason}。
      candidates: top-k 候选列表,每项 {deposit_type, confidence, evidence, source}
      reason:     空串表示正常;否则为降级原因(供前端 hint 展示)。

    流程: LLM 推理(主)→ USGS 周边矿种佐证加权 → 排序 top-k。
    USGS 只负责佐证,不直接产出矿床类型。
    """
    usgs = _usgs_nearby_context(roi_geojson)
    nearby = usgs["commodities"]

    llm_results, reason = _infer_from_llm(roi_geojson, top_k=max(top_k, 3), usgs_context=usgs)

    merged: Dict[str, Dict[str, Any]] = {}
    for name, info in llm_results.items():
        conf = info["confidence"]
        evidence = info["evidence"]
        source = info["source"]
        # 周边矿种命中 → 佐证加权
        meta = get_deposit_type_meta(name)
        commodity = meta["commodity"] if meta else ""
        if commodity and commodity in nearby:
            conf = min(1.0, conf + USGS_CORROBORATION_BONUS)
            evidence = f"{evidence};USGS MRDS 周边见同类矿种({commodity})佐证"
            source = f"{source} + USGS MRDS佐证"
        merged[name] = {
            "deposit_type": name,
            "confidence":   round(conf, 3),
            "evidence":     evidence,
            "source":       source,
        }

    sorted_list = sorted(merged.values(), key=lambda x: -x["confidence"])[:top_k]
    if not sorted_list and not reason:
        reason = "LLM 未返回有效候选,请手动选择"
    return {"candidates": sorted_list, "reason": reason}


def infer_deposit_types(roi_geojson: Dict[str, Any], top_k: int = 3) -> List[Dict[str, Any]]:
    """主入口(兼容旧调用): 仅返回 top-k 候选列表。"""
    return infer_deposit_types_detailed(roi_geojson, top_k)["candidates"]


if __name__ == "__main__":
    demos = [
        ("Chuquicamata (智利斑岩铜矿带)", {
            "type": "Polygon",
            "coordinates": [[
                [-68.95, -22.30], [-68.85, -22.30],
                [-68.85, -22.40], [-68.95, -22.40],
                [-68.95, -22.30],
            ]],
        }),
        ("胶东招远-莱州金矿带", {
            "type": "Polygon",
            "coordinates": [[
                [120.30, 37.30], [120.50, 37.30],
                [120.50, 37.45], [120.30, 37.45],
                [120.30, 37.30],
            ]],
        }),
    ]
    for title, roi in demos:
        print(f"=== {title} ===")
        res = infer_deposit_types_detailed(roi, top_k=3)
        if not res["candidates"]:
            print(f"  (降级) {res['reason']}")
        for c in res["candidates"]:
            print(f"  [{c['confidence']:.2f}] {c['deposit_type']}  ({c['source']})")
            print(f"         {c['evidence']}")
        print()
