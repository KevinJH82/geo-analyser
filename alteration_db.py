"""
蚀变-矿床知识库查询层

封装对 alteration_deposit_db.json 的所有访问:
- 矿种 / 矿床类型 两级查询(兜底手动模式用)
- 按矿床类型名查推荐蚀变矿物列表,每个矿物附带按当前传感器解析好的
  波段比值表达式 + Crosta PCA 参数
- 列出所有 deposit_type 名称(供 LLM 推理时作枚举)
"""

from __future__ import annotations

import json
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Any

_DB_PATH = Path(__file__).parent / "alteration_deposit_db.json"

# 前端传感器字符串 → JSON 中的 sensor key 映射
# 前端会先去掉 _L2 后缀,这里映射规范化后的形式
_SENSOR_MAP = {
    "Landsat8/9": "Landsat8",
    "Landsat8":   "Landsat8",
    "Sentinel2":  "Sentinel2",
    "ASTER":      "ASTER",
    "EnMAP":      "EnMAP",
    "PRISMA":     "PRISMA",
}


# ─────────────────────────────────────────────
# 油气虚拟矿种注入
# ─────────────────────────────────────────────
# JSON 里 oil_gas 是描述性结构,没有结构化 sensor_ratio。
# 这里把它转成 commodities 兼容格式,在 Python 层注入到 _load_db() 返回值里。
# 表达式依据:
#   - 红层褪色(B4/B2 取低值)→ 数学反转为 B2/B4 取高值(微渗漏区 B2/B4 偏高)
#   - 复合烃微渗漏指数:Springer 2019 (Assam-Arakan) (B2+B5)/(B3+B4) Landsat 7 → L8: (B3+B6)/(B4+B5)
#   - 粘土/绢云母:Landsat B6/B7, Sentinel-2 B11/B12, ASTER B5/B6
#   - 次生碳酸盐化:ASTER TIR B13/B14
#   - 黄铁矿化:Landsat B4/B2(铁染), ASTER B4/B5
#   - 植被红边胁迫:Sentinel-2 NDRE = (B8-B5)/(B8+B5)
#   - 热红外异常:Landsat B10 / ASTER B14 单波段亮温
#   - 高光谱烃指数 / 放射性铀:遥感不可直接算,标 unavailable

_OIL_GAS_ALTERATIONS = [
    {
        "mineral": "红层褪色异常",
        "zone": "微渗漏中心垂直投影区",
        "priority": 1,
        "anomaly_type": "Fe³⁺ → Fe²⁺ 还原褪色",
        "sensor_ratio": {
            "Landsat8":  "B2/B4",   # 反转:原 B4/B2 取低值 ≡ B2/B4 取高值
            "Sentinel2": "B2/B4",
            "ASTER":     "B1/B2",
        },
        "crosta_pca": None,
        # 红层褪色 = 烃类还原使 Fe³⁺ 减少 → 0.90µm 铁吸收"变浅",故 sign=-1(异常取浅吸收)
        "enmap_feature": {"feature_um": 0.90, "shoulder_um": [0.78, 1.05], "sign": -1},
    },
    {
        "mineral": "复合烃微渗漏指数",
        "zone": "微渗漏综合带",
        "priority": 1,
        "anomaly_type": "多波段复合(Assam-Arakan 公式)",
        "sensor_ratio": {
            "Landsat8":  "(B3+B6)/(B4+B5)",
            "Sentinel2": "(B3+B11)/(B4+B8)",
        },
        "crosta_pca": None,
    },
    {
        "mineral": "次生碳酸盐化(ΔC)",
        "zone": "微渗漏中心-边缘",
        "priority": 1,
        "anomaly_type": "CO₃²⁻ TIR 吸收",
        "sensor_ratio": {
            "ASTER": "B13/B14",
        },
        "crosta_pca": {
            "sensor": "ASTER",
            "input_bands": ["B10", "B12", "B13", "B14"],
            "pc_criterion": {
                "positive_bands": ["B13"],
                "negative_bands": ["B14"],
                "typical_pc_index": "PC2或PC3",
                "anomaly_color": "亮色=碳酸盐化",
            },
        },
        # CO₃²⁻ 在 2.30–2.35µm 的诊断吸收;EnMAP SWIR 直接测,优于 ASTER TIR 代理
        "enmap_feature": {"feature_um": 2.335, "shoulder_um": [2.26, 2.40], "sign": 1},
    },
    {
        "mineral": "粘土矿化(微渗漏诱发)",
        "zone": "微渗漏外环",
        "priority": 2,
        "anomaly_type": "Al-OH 高岭石/伊利石",
        "sensor_ratio": {
            "Landsat8":  "B6/B7",
            "Sentinel2": "B11/B12",
            "ASTER":     "B5/B6",
        },
        "crosta_pca": {
            "sensor": "ASTER",
            "input_bands": ["B1", "B4", "B6", "B7"],
            "landsat8_bands": ["B2", "B5", "B6", "B7"],
            "pc_criterion": {
                "positive_bands": ["B4"],
                "negative_bands": ["B6", "B7"],
                "typical_pc_index": "PC3或PC4",
                "anomaly_color": "亮色=粘土化",
            },
        },
        # Al-OH 在 2.20µm 的诊断吸收(高岭石/伊利石);EnMAP 窄波段精确捕捉
        "enmap_feature": {"feature_um": 2.205, "shoulder_um": [2.12, 2.245], "sign": 1},
    },
    {
        "mineral": "黄铁矿化(铁染)",
        "zone": "微渗漏带边界",
        "priority": 2,
        "anomaly_type": "Fe³⁺ 黄钾铁矾",
        "sensor_ratio": {
            "Landsat8": "B4/B2",
            "ASTER":    "B4/B5",
        },
        "crosta_pca": None,
        # Fe³⁺(黄钾铁矾/铁染)0.90µm 吸收,正异常(吸收深=富集)
        "enmap_feature": {"feature_um": 0.90, "shoulder_um": [0.78, 1.05], "sign": 1},
    },
    {
        "mineral": "植被红边胁迫(NDRE)",
        "zone": "微渗漏地表植被",
        "priority": 2,
        "anomaly_type": "烃类胁迫导致叶绿素↓",
        "sensor_ratio": {
            "Sentinel2": "(B8-B5)/(B8+B5)",   # NDRE,渗漏区取低值,这里仍用原向以 vmax/vmin 论
        },
        "crosta_pca": None,
        "note": "渗漏区 NDRE 偏低,实际算法仍按高值阈值;如需反向请人工解读",
    },
    {
        "mineral": "热红外温度异常",
        "zone": "微渗漏潜在带",
        "priority": 3,
        "anomaly_type": "LST 局部升高(轻烃挥发)",
        "sensor_ratio": {
            "Landsat8": "B10",       # 单波段亮温,异常区偏高
            "ASTER":    "B14",
        },
        "crosta_pca": None,
        "note": "单时相 LST 参考价值有限,建议结合多时相昼夜温差分析",
    },
    {
        "mineral": "高光谱烃指数(HI 1.73μm)",
        "zone": "土壤吸附烃",
        "priority": 3,
        "anomaly_type": "1.73μm CH 吸收",
        "sensor_ratio": {},   # 多光谱无对应波段
        "crosta_pca": None,
        # 土壤吸附烃在 1.73µm 的 C-H 吸收 — 只有高光谱(EnMAP)能测,多光谱无此波段
        "enmap_feature": {"feature_um": 1.73, "shoulder_um": [1.69, 1.78], "sign": 1},
        "note": "需高光谱数据(PRISMA / EnMAP / GF-5),多光谱无对应波段;EnMAP 可直接测",
    },
    {
        "mineral": "放射性铀迁移异常",
        "zone": "微渗漏顶部",
        "priority": 3,
        "anomaly_type": "U 自迁移富集",
        "sensor_ratio": {},
        "crosta_pca": None,
        "note": "需航空伽马能谱 K/U/Th 比值,不属于光学遥感",
    },
]

_OIL_GAS_INJECT = [
    {
        "commodity":    "石油",
        "commodity_en": "Oil",
        "category":     "能源矿产",
        "deposit_types": [
            {
                "type_name":         "常规油藏(微渗漏蚀变模式)",
                "type_en":           "Conventional Oil Reservoir (microseepage)",
                "tectonic_setting":  "克拉通-裂谷盆地、被动陆缘",
                "alteration_zoning": "红层褪色(中心) → 粘土化+碳酸盐化 → 黄铁矿化 → 植被胁迫(地表)",
                "host_rocks":        ["砂岩储层", "碳酸盐岩储层"],
                "ore_elements":      ["C", "H"],
                "alterations":       _OIL_GAS_ALTERATIONS,
            },
            {
                "type_name":         "致密油藏(微渗漏蚀变)",
                "type_en":           "Tight Oil (microseepage)",
                "tectonic_setting":  "克拉通盆地(鄂尔多斯/松辽)",
                "alteration_zoning": "弱化版常规油藏蚀变模式",
                "host_rocks":        ["致密砂岩", "页岩"],
                "ore_elements":      ["C", "H"],
                "alterations":       _OIL_GAS_ALTERATIONS,
            },
        ],
    },
    {
        "commodity":    "天然气",
        "commodity_en": "Natural Gas",
        "category":     "能源矿产",
        "deposit_types": [
            {
                "type_name":         "常规气藏(微渗漏蚀变模式)",
                "type_en":           "Conventional Gas Reservoir (microseepage)",
                "tectonic_setting":  "克拉通盆地、被动陆缘",
                "alteration_zoning": "同石油微渗漏,但还原性更强,黄铁矿化更明显",
                "host_rocks":        ["砂岩储层", "碳酸盐岩储层"],
                "ore_elements":      ["C", "H"],
                "alterations":       _OIL_GAS_ALTERATIONS,
            },
            {
                "type_name":         "煤层气(微渗漏蚀变)",
                "type_en":           "Coalbed Methane (microseepage)",
                "tectonic_setting":  "含煤盆地(沁水/鄂尔多斯)",
                "alteration_zoning": "蚀变信号比常规气弱,以构造解译为主",
                "host_rocks":        ["煤层", "煤系泥岩"],
                "ore_elements":      ["C", "H"],
                "alterations":       _OIL_GAS_ALTERATIONS,
            },
            {
                "type_name":         "天然气水合物(微渗漏)",
                "type_en":           "Gas Hydrate (microseepage)",
                "tectonic_setting":  "深海陆坡 / 多年冻土区",
                "alteration_zoning": "陆上多年冻土区可见 LST 异常 + 麻坑/塌陷",
                "host_rocks":        ["海底沉积", "多年冻土"],
                "ore_elements":      ["C", "H"],
                "alterations":       _OIL_GAS_ALTERATIONS,
            },
        ],
    },
]


@lru_cache(maxsize=1)
def _load_db() -> Dict[str, Any]:
    with open(_DB_PATH, "r", encoding="utf-8") as f:
        db = json.load(f)
    # 在内存中追加石油/天然气虚拟矿种(不修改源 JSON)
    db.setdefault("commodities", []).extend(_OIL_GAS_INJECT)
    return db


def normalize_sensor(sensor: str) -> Optional[str]:
    """把前端的 sensor 字符串规范化成 JSON 中使用的 key。"""
    s = sensor.replace("_L2", "")
    return _SENSOR_MAP.get(s)


# ─────────────────────────────────────────────
# 一级 / 二级查询(兜底手动模式)
# ─────────────────────────────────────────────

def list_commodities() -> List[Dict[str, str]]:
    """27 个矿种总览,用于前端一级下拉。"""
    db = _load_db()
    return [
        {
            "commodity":    c["commodity"],
            "commodity_en": c.get("commodity_en", ""),
            "category":     c.get("category", ""),
        }
        for c in db.get("commodities", [])
    ]


def list_deposit_types(commodity: str) -> List[Dict[str, str]]:
    """指定矿种下的所有矿床类型,用于前端二级下拉。"""
    db = _load_db()
    for c in db.get("commodities", []):
        if c["commodity"] == commodity:
            return [
                {
                    "type_name":         dt["type_name"],
                    "type_en":           dt.get("type_en", ""),
                    "tectonic_setting":  dt.get("tectonic_setting", ""),
                    "alteration_zoning": dt.get("alteration_zoning", ""),
                }
                for dt in c.get("deposit_types", [])
            ]
    return []


# ─────────────────────────────────────────────
# 全量 deposit_type 名称(供 LLM 枚举)
# ─────────────────────────────────────────────

@lru_cache(maxsize=1)
def all_deposit_type_names() -> List[str]:
    """所有矿种下出现过的 deposit_type 名称列表(去重保序)。"""
    db = _load_db()
    seen = set()
    out: List[str] = []
    for c in db.get("commodities", []):
        for dt in c.get("deposit_types", []):
            name = dt["type_name"]
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out


@lru_cache(maxsize=1)
def deposit_type_index() -> Dict[str, Dict[str, Any]]:
    """{deposit_type_name: deposit_type_dict_with_commodity} — 用于按矿床类型名直接查。"""
    db = _load_db()
    out: Dict[str, Dict[str, Any]] = {}
    for c in db.get("commodities", []):
        for dt in c.get("deposit_types", []):
            name = dt["type_name"]
            if name in out:
                continue
            entry = dict(dt)
            entry["_commodity"] = c["commodity"]
            entry["_commodity_en"] = c.get("commodity_en", "")
            entry["_category"] = c.get("category", "")
            out[name] = entry
    return out


# ─────────────────────────────────────────────
# 推荐蚀变目标查询(核心)
# ─────────────────────────────────────────────

# 只去除全角中文括号注释,如 "B13/B14（取低值）"
# 不能匹配半角 () —— 会误删 "(B3+B6)/(B4+B5)" 这类数学公式的括号
_ZH_PAREN = re.compile(r"（[^（）]*）")


def _clean_ratio_expr(expr: Optional[str]) -> Optional[str]:
    """清洗 ratio 表达式: 去中文括号注释、去空白。返回 None 表示不可用。"""
    if not expr:
        return None
    cleaned = _ZH_PAREN.sub("", expr).strip()
    cleaned = re.sub(r"\s+", "", cleaned)
    if not cleaned:
        return None
    # 必须是 BN 形式的表达式(不支持 PRISMA 的 R(2.10) 这类波长引用)
    if not re.search(r"B\d+", cleaned):
        return None
    return cleaned


def _build_pca_spec(crosta: Optional[Dict[str, Any]], sensor_key: str) -> Optional[Dict[str, Any]]:
    """
    根据 JSON 中的 crosta_pca 段构建该传感器下的 PCA 规格。
    返回 {input_bands, positive_bands, negative_bands, typical_pc_index, anomaly_color}
    或 None 表示不可用。
    """
    if not crosta or not isinstance(crosta, dict):
        return None

    # 选 input_bands: 默认 input_bands(通常是 ASTER 波段),如果是 Landsat8 优先用 landsat8_bands
    bands = None
    if sensor_key == "Landsat8" and crosta.get("landsat8_bands"):
        bands = crosta["landsat8_bands"]
    elif sensor_key == crosta.get("sensor"):
        bands = crosta.get("input_bands")
    elif sensor_key == "ASTER" and crosta.get("input_bands"):
        bands = crosta.get("input_bands")

    if not bands or not isinstance(bands, list):
        return None

    crit = crosta.get("pc_criterion") or {}
    pos = crit.get("positive_bands") or []
    neg = crit.get("negative_bands") or []
    if not pos and not neg:
        return None

    return {
        "input_bands":      list(bands),
        "positive_bands":   list(pos),
        "negative_bands":   list(neg),
        "typical_pc_index": crit.get("typical_pc_index", ""),
        "anomaly_color":    crit.get("anomaly_color", ""),
    }


def get_recommended_targets(deposit_type_name: str, sensor: str) -> List[Dict[str, Any]]:
    """
    输入矿床类型名 + 传感器,返回该矿床类型下所有推荐蚀变矿物的分析规格。

    每个返回项:
        {
          mineral, zone, priority, anomaly_type, absorption_um, reflectance_peak_um,
          ratio_expr,       # 已按 sensor 选定并清洗的字符串,如 "B5/B6"; None=不可用
          pca_spec,         # {input_bands, positive_bands, negative_bands, ...}; None=不可用
          ratio_available,  # bool
          pca_available,    # bool
        }

    自动过滤掉两种方法都不可用的条目。
    """
    sensor_key = normalize_sensor(sensor)
    if sensor_key is None:
        return []

    idx = deposit_type_index()
    dt = idx.get(deposit_type_name)
    if not dt:
        return []

    results = []
    for alt in dt.get("alterations", []):
        ratio_raw = (alt.get("sensor_ratio") or {}).get(sensor_key)
        ratio_expr = _clean_ratio_expr(ratio_raw)
        pca_spec = _build_pca_spec(alt.get("crosta_pca"), sensor_key)

        ratio_available = ratio_expr is not None
        pca_available = pca_spec is not None
        if not ratio_available and not pca_available:
            continue

        results.append({
            "mineral":             alt["mineral"],
            "zone":                alt.get("zone", ""),
            "priority":            int(alt.get("priority", 2)),
            "anomaly_type":        alt.get("anomaly_type", ""),
            "absorption_um":       alt.get("absorption_um"),
            "reflectance_peak_um": alt.get("reflectance_peak_um"),
            "ratio_expr":          ratio_expr,
            "pca_spec":            pca_spec,
            "ratio_available":     ratio_available,
            "pca_available":       pca_available,
        })

    # 按 priority 升序(1 在前)
    results.sort(key=lambda x: (x["priority"], x["mineral"]))
    return results


_ALL_SENSORS = ("ASTER", "Landsat8", "Sentinel2")


def get_targets_multi_sensor(deposit_type_name: str) -> List[Dict[str, Any]]:
    """
    返回矿床类型下每个推荐蚀变矿物在所有传感器上的可用情况,
    并标注首选传感器(用于自动选传感器)。

    输出每项:
      {
        mineral, zone, priority, anomaly_type,
        per_sensor: {
            "ASTER":    {ratio_expr, ratio_available, pca_spec, pca_available},
            "Landsat8": {...}, "Sentinel2": {...}
        },
        preferred_sensor: "ASTER" | "Landsat8" | "Sentinel2" | None,
      }
    """
    idx = deposit_type_index()
    dt = idx.get(deposit_type_name)
    if not dt:
        return []

    out: List[Dict[str, Any]] = []
    for alt in dt.get("alterations", []):
        per_sensor: Dict[str, Dict[str, Any]] = {}
        any_available = False
        for sk in _ALL_SENSORS:
            ratio_raw = (alt.get("sensor_ratio") or {}).get(sk)
            ratio_expr = _clean_ratio_expr(ratio_raw)
            pca_spec = _build_pca_spec(alt.get("crosta_pca"), sk)
            ra, pa = ratio_expr is not None, pca_spec is not None
            per_sensor[sk] = {
                "ratio_expr":      ratio_expr,
                "ratio_available": ra,
                "pca_spec":        pca_spec,
                "pca_available":   pa,
            }
            if ra or pa:
                any_available = True

        # 高光谱(EnMAP / PRISMA):吸收深度法,独立于多光谱 ratio/pca。
        # enmap_feature 是与传感器无关的诊断吸收特征 {feature_um, shoulder_um, sign},
        # 只要影像有对应波长即可计算,故 EnMAP / PRISMA 共用同一份特征定义。
        enmap_feature = alt.get("enmap_feature")
        bd = enmap_feature is not None
        for _hs in ("EnMAP", "PRISMA"):
            per_sensor[_hs] = {
                "ratio_expr":           None,
                "ratio_available":      False,
                "pca_spec":             None,
                "pca_available":        False,
                "band_depth_available": bd,
                "enmap_feature":        enmap_feature,
            }
        if bd:
            any_available = True

        if not any_available:
            continue

        # 首选传感器:
        #   优先 crosta_pca.sensor(JSON 注明的首选,通常 ASTER)且其 PCA 可用
        #   否则按 (ratio_available + pca_available) 总分最高,平局取 ASTER > Sentinel2 > Landsat8
        preferred = None
        crosta_sensor = (alt.get("crosta_pca") or {}).get("sensor")
        if crosta_sensor in per_sensor and per_sensor[crosta_sensor]["pca_available"]:
            preferred = crosta_sensor

        if preferred is None:
            ranked = sorted(
                _ALL_SENSORS,
                key=lambda s: -(int(per_sensor[s]["ratio_available"]) + int(per_sensor[s]["pca_available"])),
            )
            for s in ranked:
                if per_sensor[s]["ratio_available"] or per_sensor[s]["pca_available"]:
                    preferred = s
                    break

        out.append({
            "mineral":             alt["mineral"],
            "zone":                alt.get("zone", ""),
            "priority":            int(alt.get("priority", 2)),
            "anomaly_type":        alt.get("anomaly_type", ""),
            "absorption_um":       alt.get("absorption_um"),
            "reflectance_peak_um": alt.get("reflectance_peak_um"),
            "enmap_feature":       alt.get("enmap_feature"),
            "per_sensor":          per_sensor,
            "preferred_sensor":    preferred,
        })

    out.sort(key=lambda x: (x["priority"], x["mineral"]))
    return out


def get_deposit_type_meta(deposit_type_name: str) -> Optional[Dict[str, Any]]:
    """返回矿床类型的元信息(矿种归属、大地构造背景、蚀变分带描述)。"""
    dt = deposit_type_index().get(deposit_type_name)
    if not dt:
        return None
    return {
        "type_name":         dt["type_name"],
        "type_en":           dt.get("type_en", ""),
        "commodity":         dt.get("_commodity", ""),
        "commodity_en":      dt.get("_commodity_en", ""),
        "category":          dt.get("_category", ""),
        "tectonic_setting":  dt.get("tectonic_setting", ""),
        "alteration_zoning": dt.get("alteration_zoning", ""),
        "ore_elements":      dt.get("ore_elements", []),
        "host_rocks":        dt.get("host_rocks", []),
    }


if __name__ == "__main__":
    print(f"矿种数: {len(list_commodities())}")
    print(f"所有 deposit_type 数: {len(all_deposit_type_names())}")
    print()
    targets = get_recommended_targets("斑岩型铜矿", "Landsat8/9")
    print(f"斑岩型铜矿 / Landsat8 推荐蚀变目标({len(targets)} 个):")
    for t in targets:
        ra = "✓" if t["ratio_available"] else "✗"
        pa = "✓" if t["pca_available"] else "✗"
        print(f"  [P{t['priority']}] {t['mineral']:<12} ratio={ra}({t['ratio_expr']}) pca={pa}")
