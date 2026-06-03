"""
insar_broker.py — 订阅 geo-insar 标准输出目录(Phase 1.5)

geo-analyser 与 geo-insar 都在同机部署,通过文件系统订阅最简单:
- 扫描 /opt/deepexplor-services/geo-insar/downloads/ 下的所有 AOI 目录
- 每个 AOI 若有 sentinel1_insar/<pair>/metadata.json,认为是可分析的堆栈
- 提供 list_available_stacks() 给前端列表,以及 get_stack_path() 给分析模块用

不引入额外的消息队列,纯文件系统订阅,失败容忍度高。
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional

# 与 commons/insar_schema.json 对齐的默认路径
DEFAULT_GEO_INSAR_OUTPUTS = "/opt/deepexplor-services/geo-insar/downloads"


def scan_available_aois(geo_insar_outputs: str = DEFAULT_GEO_INSAR_OUTPUTS) -> List[Dict]:
    """
    扫描 geo-insar 输出目录,返回所有可分析的 AOI 列表。

    Returns
    -------
    [
        {
            "aoi_name": str,
            "aoi_path": str,
            "n_pairs": int,
            "date_range": (earliest, latest),
            "stack_index_path": str | None,
        },
        ...
    ]
    """
    root = Path(geo_insar_outputs)
    if not root.exists():
        return []

    out = []
    for aoi_dir in root.iterdir():
        if not aoi_dir.is_dir():
            continue
        pair_dirs = sorted(aoi_dir.glob("sentinel1_insar/*"))
        pair_dirs = [p for p in pair_dirs if p.is_dir() and (p / "metadata.json").exists()]
        if not pair_dirs:
            continue

        # 优先用 stack_index.json(汇总),fallback 到逐个 pair
        idx = aoi_dir / "stack_index.json"
        if idx.exists():
            try:
                with open(idx, "r", encoding="utf-8") as f:
                    summary = json.load(f)
                out.append({
                    "aoi_name": aoi_dir.name,
                    "aoi_path": str(aoi_dir),
                    "n_pairs": summary.get("pair_count", len(pair_dirs)),
                    "date_range": summary.get("date_range", [None, None]),
                    "stack_index_path": str(idx),
                })
                continue
            except Exception:
                pass

        # fallback
        dates = []
        for pdir in pair_dirs:
            try:
                with open(pdir / "metadata.json", "r", encoding="utf-8") as f:
                    m = json.load(f)
                dates.extend([m.get("master_date"), m.get("slave_date")])
            except Exception:
                pass
        dates = sorted(d for d in dates if d)
        out.append({
            "aoi_name": aoi_dir.name,
            "aoi_path": str(aoi_dir),
            "n_pairs": len(pair_dirs),
            "date_range": [dates[0], dates[-1]] if dates else [None, None],
            "stack_index_path": None,
        })
    return out


def get_stack_path(aoi_name: str,
                   geo_insar_outputs: str = DEFAULT_GEO_INSAR_OUTPUTS) -> Optional[str]:
    """根据 aoi_name 返回对应 AOI 目录路径(供 insar_timeseries.load_insar_stack 用)。"""
    aoi_dir = Path(geo_insar_outputs) / aoi_name
    if aoi_dir.is_dir():
        return str(aoi_dir)
    return None
