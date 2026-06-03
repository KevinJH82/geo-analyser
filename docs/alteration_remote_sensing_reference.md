# 多矿种 × 多卫星 蚀变遥感对照手册

> 用途:对照核校 `alteration_deposit_db.json` 中各矿种的蚀变矿物 / 波段比值表达式 / Crósta PCA 输入波段;为扩充新矿种或新传感器提供文献依据。
> 
> 说明:以下所有公式整理自下列经典文献,具体使用前建议核对原文细节(不同作者对同名指数的版本可能有差异)。我已为每条公式标注最高频引用来源,但**任何投入生产的指数表达式都应交叉验证后再入库**。

---

## 一、卫星传感器波段速查

### ASTER(VNIR 3 + SWIR 6 + TIR 5 = 14 波段)
| 波段 | 范围(μm) | 分辨率 | 用途 |
|---|---|---|---|
| B1 | 0.52–0.60 | 15 m | Green VNIR |
| B2 | 0.63–0.69 | 15 m | Red VNIR |
| B3N / B3B | 0.78–0.86 | 15 m | NIR(Nadir / Backward,本系统统一用 B3N → B3) |
| B4 | 1.60–1.70 | 30 m | SWIR1 |
| B5 | 2.145–2.185 | 30 m | SWIR2.16 |
| B6 | 2.185–2.225 | 30 m | SWIR2.20(Al-OH 强吸收) |
| B7 | 2.235–2.285 | 30 m | SWIR2.26 |
| B8 | 2.295–2.365 | 30 m | SWIR2.33(CO3 强吸收) |
| B9 | 2.36–2.43 | 30 m | SWIR2.40 |
| B10–14 | 8.125–11.65 | 90 m | TIR(石英、长石、碳酸盐 TIR 特征) |

### Sentinel-2(MSI,13 波段)
| 波段 | 中心(μm) | 分辨率 | 用途 |
|---|---|---|---|
| B2 | 0.490 | 10 m | Blue |
| B3 | 0.560 | 10 m | Green |
| B4 | 0.665 | 10 m | Red |
| B5 | 0.705 | 20 m | Red Edge 1 |
| B6 | 0.740 | 20 m | Red Edge 2 |
| B7 | 0.783 | 20 m | Red Edge 3 |
| B8 | 0.842 | 10 m | NIR(宽) |
| B8A | 0.865 | 20 m | NIR(窄) |
| B11 | 1.610 | 20 m | SWIR1 |
| B12 | 2.190 | 20 m | SWIR2(覆盖 Al-OH 2.20μm,缺 CO3 2.33μm) |

> **关键限制**:Sentinel-2 SWIR 只有 B11/B12 两个,无法区分黏土族细分(无法分辨绢云母 vs 高岭石 vs 明矾石),也不能直接做碳酸盐 2.33μm 吸收检测。

### Landsat 8 / 9(OLI 9 波段)
| 波段 | 中心(μm) | 分辨率 | 用途 |
|---|---|---|---|
| B1 | 0.443 | 30 m | Coastal Aerosol |
| B2 | 0.482 | 30 m | Blue |
| B3 | 0.561 | 30 m | Green |
| B4 | 0.655 | 30 m | Red |
| B5 | 0.865 | 30 m | NIR |
| B6 | 1.609 | 30 m | SWIR1 |
| B7 | 2.201 | 30 m | SWIR2(覆盖 Al-OH,但带宽 ~180 nm,与 ASTER B5–9 相比无法细分黏土族) |
| B10–11 | 10.9 / 12.0 | 100 m | TIR |

---

## 二、常见蚀变矿物的光谱诊断

| 矿物 | 中心吸收(μm) | 元素 | 备注 |
|---|---|---|---|
| 赤铁矿 (Hematite) | 0.55 / 0.85–0.87 | Fe³⁺ | 红层、铁帽 |
| 褐铁矿 (Goethite) | 0.48 / 0.90 | Fe³⁺ | 氧化带 |
| 黄钾铁矾 (Jarosite) | 0.43 / 0.93 | Fe³⁺, K, SO₄ | 高级泥化指示 |
| 绢云母 / 白云母 | 2.20(Al-OH) | K, Al | 斑岩 phyllic 带核心 |
| 伊利石 | 2.20(Al-OH) | K, Al | 与绢云母光谱难分 |
| 蒙脱石 | 2.20(Al-OH) + 1.91(H₂O) | Al, Si | argillic 带 |
| 高岭石 | 2.16 + 2.20 双峰(Al-OH) | Al | argillic 带 |
| 明矾石 | 2.17(Al-OH) + 1.76(SO₄) | K, Al, S | advanced argillic 带 |
| 叶蜡石 | 2.165(Al-OH) | Al | advanced argillic 高温端 |
| 绿泥石 | 2.25(Fe-OH) + 2.33(Mg-OH) | Fe, Mg | propylitic 带 |
| 绿帘石 | 2.34(Mg-OH) | Ca, Al, Fe | propylitic 带 |
| 方解石 / 白云石 | 2.33–2.35(CO₃) | Ca, Mg | 碳酸盐化、矽卡岩外带 |
| 滑石 | 2.32(Mg-OH) | Mg | 超基性蚀变 |
| 蛇纹石 | 2.32(Mg-OH) | Mg | 超基性蚀变 |
| 石膏 | 1.45 / 1.94 / 2.21 | Ca, SO₄ | 蒸发岩 / 蚀变伴生 |

---

## 三、按矿种 × 卫星 × 方法 对照表

> 表格说明:
> - **波段比值法** 列给出代数表达式,数值越高 = 该蚀变信号越强(已校核方向)
> - **Crósta PCA** 列给出推荐输入波段子集,目标 PC 由载荷符号自动选(本系统 `alteration_analysis.calc_crosta_pca` 已实现自动符号判定,等同 Loughlin 1991 改进版)
> - 仅当卫星具备对应吸收波段才支持该矿物;无可用方法的格子标 "—"

### 3.1 斑岩型铜矿(Porphyry Cu / Cu-Au / Cu-Mo)

蚀变 zoning(由内向外):钾化(黑云母-钾长石)→ 绢云母化(phyllic)→ 泥化(argillic)→ 青磐岩化(propylitic)

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | Landsat 8 PCA 输入 | Sentinel-2 比值 | 主要文献 |
|---|---|---|---|---|---|---|---|---|
| 绢云母化(phyllic,绢云母 / 伊利石) | 1 | K, Al | `B4/B6` 或 `(B5+B7)/B6` | B1, B4, B6, B7 | `B6/B7` | B2, B5, B6, B7 | `B11/B12`(粗略) | Sabins 1999;Mars & Rowan 2006 |
| 高岭石(argillic) | 1 | Al | `B4/B5` 或 `(B4+B6)/B5` | B1, B4, B5, B6 | `B6/B7`(无法细分) | B2, B5, B6, B7 | `B11/B12` | Mars & Rowan 2006 |
| 明矾石(advanced argillic) | 2 | K, Al, S | `B7/B6` | B1, B3, B6, B7 | — | — | — | Crósta et al. 2003 |
| 叶蜡石 | 2 | Al | `B5/B6`(2.165 吸收) | B1, B5, B6, B7 | — | — | — | Hunt 1977;Mars 2018 |
| 绿泥石(propylitic) | 2 | Fe, Mg | `B5/B8` 或 `(B6+B9)/(B7+B8)` | B1, B5, B7, B8 | `B7/B5`(粗) | B2, B5, B6, B7 | `B12/B8A` | Rowan & Mars 2003 |
| 绿帘石(propylitic) | 2 | Ca, Fe, Al | `B7/B8` 或 `(B6+B9)/(B7+B8)` | B1, B6, B8, B9 | — | — | — | Rowan & Mars 2003 |
| 铁帽 / 黄铁矿氧化(赤铁矿) | 1 | Fe³⁺ | `B2/B1` | B1, B2, B3, B4 | `B4/B2` | B2, B4, B5, B6 | `B4/B2` | Sabins 1999 |

### 3.2 浅成低温热液型金矿(Epithermal Au / Au-Ag)

高硫型蚀变 zoning:硅化核 → 高级泥化(明矾石-叶蜡石-高岭石) → 绢云母化 → 青磐岩化

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | Landsat 8 PCA 输入 | Sentinel-2 比值 | 主要文献 |
|---|---|---|---|---|---|---|---|---|
| 明矾石 | 1 | K, Al, S | `B7/B6` | B1, B3, B6, B7 | — | — | — | Crósta et al. 2003;Bedini 2009 |
| 叶蜡石 | 1 | Al | `B5/B6` | B1, B5, B6, B7 | — | — | — | Mars & Rowan 2006 |
| 高岭石 / 地开石 | 1 | Al | `B4/B5` | B1, B4, B5, B6 | `B6/B7` | B2, B5, B6, B7 | `B11/B12` | Cudahy 1997 |
| 绢云母 / 伊利石 | 2 | K, Al | `B4/B6` 或 `(B5+B7)/B6` | B1, B4, B6, B7 | `B6/B7` | B2, B5, B6, B7 | `B11/B12` | Sabins 1999 |
| 硅化 / 蛋白石 | 2 | Si | `B14/B12`(TIR Reststrahlen) | TIR 单做 | — | — | — | Rowan & Mars 2003 |
| 黄铁矿氧化(jarosite) | 1 | Fe, K, SO₄ | `B2/B1` 配 `B5/B4` | B1, B2, B3, B4 | `B4/B2` 配 `B7/B5` | B2, B4, B5, B6 | `B4/B2` | Loughlin 1991 |

### 3.3 造山型金矿(Orogenic Au,绿岩带 / 变质流体型)

主要蚀变:绢云母化 + 碳酸盐化(CO₂ 流体)+ 绿泥石化 + 黄铁矿化

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | Landsat 8 PCA 输入 | Sentinel-2 比值 | 主要文献 |
|---|---|---|---|---|---|---|---|---|
| 绢云母化 | 1 | K, Al | `B4/B6` | B1, B4, B6, B7 | `B6/B7` | B2, B5, B6, B7 | `B11/B12` | Yousefi & Carranza 2015 |
| 铁白云石 / 菱铁矿(碳酸盐化) | 1 | Ca, Mg, Fe | `B8/B6` 或 `(B6+B9)/(B7+B8)` | B1, B3, B7, B8 | — | — | — | Rockwell 2013 |
| 绿泥石化 | 2 | Fe, Mg | `B5/B8` | B1, B5, B7, B8 | `B7/B5` | B2, B5, B6, B7 | `B12/B8A` | Crósta 1989 |
| 黄铁矿氧化(赤铁矿/褐铁矿) | 2 | Fe³⁺ | `B2/B1` | B1, B2, B3, B4 | `B4/B2` | B2, B4, B5, B6 | `B4/B2` | Sabins 1999 |

### 3.4 矽卡岩型铁/钨/铜矿(Skarn Fe / W / Cu)

矽卡岩带通常呈带状(石榴石-辉石内带,符山石-绿帘石外带),叠加大理岩化 / 角岩化

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | Sentinel-2 比值 | 主要文献 |
|---|---|---|---|---|---|---|---|
| 绿帘石 / 绿泥石(外带) | 1 | Ca, Fe, Mg | `B7/B8` 或 `(B6+B9)/(B7+B8)` | B1, B6, B8, B9 | `B7/B5` | `B12/B8A` | Mars 2018 |
| 大理岩 / 方解石(围岩) | 1 | Ca | `B8/B6` 或 `B7/B8` | B1, B3, B7, B8 | — | — | Rowan & Mars 2003 |
| 石榴石-辉石(内带,Fe-Mg 硅酸盐) | 2 | Fe, Ca, Mg | `B5/B4`(铁信号) + TIR | B1, B4, B5, B6 | `B5/B4` | `B8/B4` | van der Meer 2012 |
| 赤铁矿(矽卡岩 Fe 矿核心) | 1 | Fe³⁺ | `B2/B1` | B1, B2, B3, B4 | `B4/B2` | `B4/B2` | Sabins 1999 |

### 3.5 块状硫化物矿(VMS / SEDEX 铅锌)

蚀变:硅化-绢云母化筒(下盘) + 绿泥石化(海底蚀变)+ 黄铁矿氧化形成铁帽

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | 主要文献 |
|---|---|---|---|---|---|---|
| 绢云母化 / 伊利石 | 1 | K, Al | `B4/B6` | B1, B4, B6, B7 | `B6/B7` | Pour & Hashim 2014 |
| 绿泥石化(下盘) | 1 | Fe, Mg | `B5/B8` | B1, B5, B7, B8 | `B7/B5` | Pour & Hashim 2015 |
| 铁帽(gossan,jarosite + goethite + hematite) | 1 | Fe³⁺, K, SO₄ | `B2/B1` + `B5/B4` | B1, B2, B3, B4 | `B4/B2` + `B7/B5` | Loughlin 1991 |

### 3.6 红土型镍矿(Lateritic Ni / Ni-Co)

风化壳:腐岩带(铁帽,赤铁矿/褐铁矿) → 蛇纹石残留带 → 含镍滑石 / 蛇纹石 → 母岩

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | 主要文献 |
|---|---|---|---|---|---|
| 铁帽(赤铁矿 + 褐铁矿) | 1 | Fe³⁺, Ni | `B2/B1` + `B5/B4` | B1, B2, B3, B4 | Crósta et al. 2003 |
| 蛇纹石 / 滑石(2.32 Mg-OH) | 1 | Mg, Ni | `B7/B8`(注意 Mg-OH 在 ASTER B8 附近) | B1, B6, B8, B9 | Rowan & Mars 2003 |
| 高岭石(底部) | 2 | Al | `B4/B5` | B1, B4, B5, B6 | Mars 2018 |

### 3.7 IOCG(铁氧化物-铜-金型)

蚀变:钠化(钠长石化) → 钾化(钾长石 + 黑云母) → 钙硅酸盐(绿帘石 + 绿泥石) → 赤铁矿化(核心)

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | Landsat 8 比值 | 主要文献 |
|---|---|---|---|---|---|---|
| 赤铁矿化(核心) | 1 | Fe³⁺ | `B2/B1` | B1, B2, B3, B4 | `B4/B2` | Skirrow 2010 |
| 钾化(黑云母-钾长石) | 1 | K | `B4/B6` + TIR | B1, B4, B6, B7 | `B6/B7` | Rowan 2006 |
| 钠长石化(albitization) | 2 | Na, Al | TIR `B14/B13`(长石族) | — | — | Rowan 2006 |
| 绿帘石 / 绿泥石 | 2 | Ca, Fe, Mg | `B7/B8` | B1, B6, B8, B9 | `B7/B5` | Mars 2018 |

### 3.8 油气微渗漏蚀变(常规油藏 / 致密油气)

地表蚀变:红层褪色 + 次生碳酸盐化 + 黏土矿化 + 植被红边胁迫 + 高光谱 HI 指数

| 蚀变 | 优先级 | 元素 / 指示 | Landsat 8 比值(主推) | ASTER 比值 | Sentinel-2 比值 | 主要文献 |
|---|---|---|---|---|---|---|
| 红层褪色(铁还原为 Fe²⁺,显色减弱) | 1 | Fe²⁺/Fe³⁺ 比 | `B2/B4`(蓝/红反指数) | `B1/B2` | `B2/B4` | Saunders et al. 1999 |
| 复合烃微渗漏指数 | 1 | 综合指示 | `(B3+B6)/(B4+B5)` | — | 不直接对应 | Saunders 1999 |
| 次生碳酸盐化(地表碳酸盐胶结) | 1 | CaCO₃ | — | `B8/B6` 或 `B7/B8` | — | Yang et al. 2000 |
| 黏土矿化(伊利石 / 高岭石) | 2 | Al | `B6/B7` | `B4/B6` | `B11/B12` | Khan & Jacobson 2008 |
| 黄铁矿化(后期氧化) | 3 | Fe, S | `B4/B2` | `B2/B1` + `B5/B4` | `B4/B2` | Saunders 1999 |
| 植被红边胁迫 NDRE | 1 | 叶绿素衰减 | — | — | `(B8A−B5)/(B8A+B5)` | Noomen et al. 2008 |
| 高光谱烃指数 HI | 2 | C-H 1.73μm | — | 通常需高光谱(HyMap/PRISMA) | — | Cloutis 1989 |
| 热红外异常 | 3 | 地表温度 | `B10`(LST) | TIR B10–14 | — | Saunders 1999 |
| 放射性铀(土壤 ²²²Rn 富集副产) | 3 | U, Rn | 需航磁/航放支持 | — | — | Pirajno 2009 |

> **油气蚀变要点**:红层褪色和次生碳酸盐化是公认的两大近地表指示,但单一指数易受地表干扰(植被、土壤湿度、人工地物)误判;实际项目中应**多指数复合 + 已知井位标定**。本系统 `alteration_deposit_db.json` 的"石油 / 天然气"两个 commodity 即按此模式组织。

### 3.9 砂岩型铀矿 / 不整合面型铀矿

蚀变:层间氧化(铁氧化物 + 褐铁矿)+ 黏土化(高岭石 / 蒙脱石)+ 沥青铀矿(还原带)

| 蚀变带 / 矿物 | 优先级 | 元素 | ASTER 比值 | Landsat 8 比值 | 主要文献 |
|---|---|---|---|---|---|
| 层间氧化(赤铁矿) | 1 | Fe³⁺ | `B2/B1` | `B4/B2` | Sabins 1999 |
| 黏土化(高岭石/蒙脱石) | 1 | Al | `B4/B5` 或 `B4/B6` | `B6/B7` | Crósta et al. 2003 |
| 褪色带(还原界面) | 1 | Fe²⁺ | `B1/B2` | `B2/B4` | Saunders 1999(借鉴) |
| 沥青铀矿(放射性) | 3 | U | 需航放 | — | — |

### 3.10 其它矿种紧凑表(对照本系统数据库)

> 下列矿种 / 矿床类型已在 `alteration_deposit_db.json` 中入库,本节按"矿种 → 矿床类型 → 关键蚀变"紧凑组织,公式与数据库保持一致。需要细分文献溯源时,可逐条与 3.1–3.9 节的"经典版本"对照(参见本文档末尾的"算法差异清单")。

#### 3.10.1 铅锌(Pb-Zn)
| 矿床类型 | 蚀变 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | 主要文献 |
|---|---|---|---|---|---|---|
| MVT 型铅锌 | 白云石 | 1 | CO₃ | `B13/B14`(TIR) | B10, B12, B13, B14 | Leach et al. 2010 |
| MVT 型铅锌 | 硅化 | 2 | Si | TIR `B13` 或 `B14` 反射差 | B10, B12, B13, B14 | Hunt & Salisbury 1976 |
| MVT 型铅锌 | 黏土化 | 3 | Al | `B4/B6` | B1, B4, B6, B7 | 与斑岩 Cu 黏土公式一致 |
| SEDEX 型铅锌 | 电气石(B-OH) | 1 | B, OH | `B6/B5` | B1, B4, B5, B6 | Slack 1996 |
| SEDEX 型铅锌 | 绢云母 + 硅化 + 碳酸盐 | 1-3 | K/Al/Si/CO₃ | 同斑岩 Cu / 矽卡岩 | — | — |
| 矽卡岩型铅锌(锰质矽卡岩) | 蔷薇辉石 / 锰钙铁辉石 | 1 | Mn(TIR) | TIR 难以单波段提取 | B10, B12, B13, B14(pos=B10, neg=B14) | Cui et al. 2019 |

#### 3.10.2 金 — 补充类型(浅低硫 / 卡林型 / 斑岩 Cu-Au)
| 矿床类型 | 蚀变 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | 主要文献 |
|---|---|---|---|---|---|---|
| 浅低温热液(低硫型) | 冰长石(Adularia) | 1 | K(TIR) | TIR Reststrahlen | B10, B12, B13, B14(pos=B13, neg=B12) | Hauff 2008 |
| 浅低温热液(低硫型) | 玉髓状硅化 | 1 | Si | — | B10, B12, B13, B14(pos=B12, neg=B14) | Rowan & Mars 2003 |
| 卡林型金矿 | 似碧玉岩硅化 | 1 | Si | — | TIR PCA | Hofstra & Cline 2000 |
| 卡林型金矿 | 伊利石 + 高岭石 | 1 | Al | `B5/B6` | B1, B4, B6, B7 | Mars 2018 |
| 卡林型金矿 | 方解石负异常(去碳酸盐化) | 1 | CO₃ 缺失 | `B13/B14` 取**低**值 | B10, B12, B13, B14(pos=B13/B14, neg=B12) | Rowan & Bowers 1995 |
| 斑岩型铜金 | 同斑岩 Cu(绢云母 + 绿泥石 + 绿帘石)+ 磁铁矿氧化(`B4/B2` L8) | 1-2 | K/Fe/Mg/Fe³⁺ | 见 3.1 | — | Sillitoe 2010 |

#### 3.10.3 银(Ag)— 浅成低温热液型银矿
| 蚀变 | 优先级 | 元素 | ASTER 比值 | ASTER PCA 输入 | 备注 |
|---|---|---|---|---|---|
| 绢云母 + 冰长石 + 硅化 | 1 | K, Al, Si | 同浅低温低硫 Au | 同 3.10.2 | 蚀变模式与浅低温 Au 几乎一致,Ag 富集在 phyllic-adularia 核 |
| 锰氧化物(外环) | 2 | Mn(TIR) | — | B10, B12, B13, B14(pos=B10, neg=B14) | Hewett 1968 |

#### 3.10.4 铂族(PGE)— 层状镁铁质侵入体
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 蛇纹石 | 1 | `B6/B8` | B1, B4, B8, B9 | Rajendran & Nasir 2014 |
| 滑石 | 2 | `B6/B8` | B1, B4, B8, B9 | 同上 |
| 次闪石(Cummingtonite) | 2 | `B7/B8` | B1, B4, B8, B9 | — |

#### 3.10.5 铁(Fe)— 补充类型
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | Landsat 8 比值 | Sentinel-2 比值 | 文献 |
|---|---|---|---|---|---|---|
| BIF(条带状铁建造) | 赤铁矿/磁铁矿(→褐铁矿) | 1 | `B3/B1` | `B4/B2` | `B4/B2` | Sabins 1999 |
| BIF | 硅化条带 | 2 | TIR | — | — | — |
| 矽卡岩型铁 | 石榴子石 / 透辉石(内带) | 1 | TIR | — | — | Mars 2018 |
| 矽卡岩型铁 | 阳起石 / 绿帘石(外带) | 2 | `B7/B8` | — | — | Rowan & Mars 2003 |

#### 3.10.6 钼(Mo)
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|---|
| 斑岩型钼(Climax 型) | 钾长石(钾化强) | 1 | TIR | B10, B12, B13, B14(pos=B13, neg=B12) | Carten et al. 1993 |
| 斑岩型钼 | 硅化(强) | 1 | TIR | TIR | — |
| 斑岩型钼 | 绢云母化 | 2 | `B5/B6` | B1, B4, B6, B7 | — |
| 斑岩型铜钼 | 同斑岩 Cu,钾化带 Mo 富集更高 | 1-2 | 见 3.1 | — | Sillitoe 2010 |

#### 3.10.7 钨(W)/ 锡(Sn)
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|---|
| 矽卡岩型钨 | 石榴子石(TIR) | 1 | TIR | B10, B12, B13, B14(pos=B12, neg=B13) | Meinert et al. 2005 |
| 矽卡岩型钨 | 白云母 | 1 | `B5/B6` | B1, B4, B6, B7 | — |
| 矽卡岩型钨 | 萤石 | 2 | TIR | B10, B12, B13, B14(pos=B10, neg=B12) | Hunt & Salisbury 1976 |
| 云英岩型 W-Sn | 白云母 + 硅化 + 萤石 + 电气石 | 1-2 | 见各矿物 | 同上 | Černý et al. 2005 |

#### 3.10.8 镍(Ni)— 岩浆硫化物型(补充红土型见 3.6)
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 蛇纹石 | 1 | `B6/B8` | B1, B4, B8, B9 | Rajendran & Nasir 2014 |
| 滑石 | 2 | `B6/B8` | B1, B4, B8, B9 | — |
| 绿泥石(海底蚀变伴生) | 3 | `B6/B8` | B1, B4, B8, B9 | — |

#### 3.10.9 钴(Co)
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | 文献 |
|---|---|---|---|---|
| 沉积岩容矿铜钴 | 硅化 + 碳酸盐 + 绢云母 | 1-3 | 同斑岩 Cu | Hitzman et al. 2005 |
| 红土型钴 | 含钴针铁矿(Fe³⁺) | 1 | L8 `B4/B2` | Berger 2011 |
| 红土型钴 | 锰氧化物 | 2 | TIR PCA(pos=B10, neg=B14) | Hewett 1968 |

#### 3.10.10 锰(Mn)— 火山沉积型 / 沉积型
| 蚀变 | 优先级 | ASTER 比值 | 文献 |
|---|---|---|---|
| 硅化(围岩) | 2 | TIR | — |
| 绿泥石(海相伴生) | 3 | `B6/B8` | — |

> 沉积型锰矿主要表现为锰氧化物层位,遥感主要靠地表 Mn-Fe 氧化壳;直接可识别性弱,需结合地球化学。

#### 3.10.11 铬(Cr)— 豆荚状铬铁矿(蛇绿岩型)
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 蛇纹石(围岩蚀变) | 1 | `B6/B8` | B1, B4, B8, B9 | Rajendran et al. 2011 |

#### 3.10.12 稀土(REE)
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|---|
| 碳酸岩型稀土 | 碳酸盐(方解石/白云石) | 1 | `B13/B14`(TIR) | B10, B12, B13, B14(pos=B12, neg=B13/B14) | Bedini 2009 |
| 碳酸岩型稀土 | 赤铁矿化 | 1 | L8 `B4/B2` | — | Verplanck et al. 2014 |
| 碳酸岩型稀土 | 霓石 / 钠辉石 / 钠闪石(碱性指示) | 2 | TIR | B10, B12, B13, B14(pos=B10, neg=B12) | Mars & Rowan 2011 |
| 离子吸附型稀土 | 高岭石 + 埃洛石 | 1 | `B4/B6` | B1, B4, B6, B7 | Bao & Zhao 2008 |
| 离子吸附型稀土 | 褐铁矿 | 3 | L8 `B4/B2` | — | — |

#### 3.10.13 锂(Li)
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|---|
| 伟晶岩型(锂辉石) | 白云母(围岩) | 1 | `B5/B6` | B1, B4, B6, B7 | Cardoso-Fernandes et al. 2019 |
| 伟晶岩型 | 电气石(B-OH) | 2 | `B6/B5` | B1, B4, B5, B6 | Slack 1996 |
| 伟晶岩型 | 钠长石化 | 2 | TIR | B10, B12, B13, B14(pos=B12, neg=B13) | — |
| 伟晶岩型 | 锂辉石(2.20/1.91/2.38 μm 三重峰) | 1 | **需高光谱**(PRISMA / EnMAP) | — | Cardoso-Fernandes et al. 2021 |
| 沉积型(粘土型 Li) | 蒙脱石 / 伊利石(含锂) | 1 | `B4/B6` | B1, B4, B6, B7 | Kesler et al. 2012 |

#### 3.10.14 铀(U)— 补充类型
| 矿床类型 | 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|---|
| 不整合面型铀 | 绿泥石 + 伊利石 + 硅化 + 赤铁矿(四者共存为关键判据) | 1 | 各见 3.1–3.4 | 各见单矿物 | Jefferson et al. 2007 |
| 钠交代型铀 | 钠长石化(Na 富集) | 1 | TIR | B10, B12, B13, B14(pos=B10, neg=B12) | Cuney et al. 2012 |
| 钠交代型铀 | 绿泥石化 | 2 | `B6/B8` | B1, B4, B8, B9 | — |

#### 3.10.15 铝土矿(Bauxite)— 红土型
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 三水铝石 / 勃姆石(2.20 + 2.30 μm) | 1 | `B4/B6` | B1, B4, B6, B7 | Bárdossy 1982 |
| 高岭石(残积) | 2 | `B4/B6` | B1, B4, B6, B7 | — |
| 针铁矿 / 赤铁矿(铁帽) | 1 | L8 `B4/B2` | — | Sabins 1999 |

#### 3.10.16 锑(Sb)/ 汞(Hg)— 低温热液型
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 硅化(玉髓 / 石英) | 1 | TIR | B10, B12, B13, B14(pos=B12, neg=B14) | Hunt & Salisbury 1976 |
| 碳酸盐(方解石化) | 2 | `B13/B14` | B10, B12, B13, B14(pos=B12, neg=B13/B14) | — |
| 粘土化(高岭石) | 3 | `B4/B6` | B1, B4, B6, B7 | — |

#### 3.10.17 金刚石(Diamond)— 金伯利岩型
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 蛇纹石化(金伯利岩岩管核心) | 1 | `B6/B8` | B1, B4, B8, B9 | Mitchell 1986 |
| 碳酸盐(后期热液) | 1 | `B13/B14` | B10, B12, B13, B14 | — |
| 滑石 | 2 | `B6/B8` | B1, B4, B8, B9 | — |

> 金伯利岩管直径仅几十-几百米,需高空间分辨率影像(ASTER 15 m VNIR 仍勉强);常以**圆形 / 环形蚀变晕**为判别特征。

#### 3.10.18 萤石(Fluorite)— 热液脉型
| 蚀变 | 优先级 | ASTER 比值 | ASTER PCA 输入 | 文献 |
|---|---|---|---|---|
| 硅化(围岩) | 1 | TIR | B10, B12, B13, B14(pos=B12, neg=B14) | — |
| 萤石(TIR Reststrahlen) | 1 | TIR | B10, B12, B13, B14(pos=B10, neg=B12) | Hunt & Salisbury 1976 |
| 绢云母(伴生) | 3 | `B5/B6` | B1, B4, B6, B7 | — |

---

## 四、Crósta PCA 标准输入波段速查(经典组合)

不论矿种,Crósta PCA 的"标准 4 波段子集"原则:**1 个反射高的对照波段 + 2-3 个目标吸收带覆盖波段**,让目标矿物的吸收特征落到某一个 PC 上独立呈现。

### ASTER
| 目标 | 推荐输入 | 期望 PC 载荷符号 |
|---|---|---|
| 铁氧化物(Fe³⁺) | B1, B2, B3, B4 | B1 负 / B2 正 |
| Al-OH 黏土(广义) | B1, B4, B6, B7 | B4 正 / B6 负 |
| 高岭石 vs 绢云母细分 | B1, B4, B5, B6 / B1, B4, B6, B7 | 同上,选峰位 |
| 明矾石 | B1, B3, B6, B7 | B6 负 / B7 正 |
| 碳酸盐(CO₃) | B1, B3, B7, B8 | B7 正 / B8 负 |
| Mg-OH(绿泥石 / 绿帘石 / 滑石) | B1, B5, B7, B8 | B5 正 / B8 负 |

### Landsat 8(注:OLI 缺细分能力,主要靠铁 / 黏土二分)
| 目标 | 推荐输入 | 期望 PC 载荷符号 |
|---|---|---|
| 铁氧化物 | B2, B4, B5, B6 | B2 负 / B4 正 |
| 含羟基矿物(广义黏土) | B2, B5, B6, B7 | B6 正 / B7 负 |

### Sentinel-2
| 目标 | 推荐输入 | 备注 |
|---|---|---|
| 铁氧化物 | B2, B3, B4, B8 | B2 负 / B4 正 |
| 含羟基矿物(粗略) | B2, B4, B11, B12 | B11 正 / B12 负 |

---

## 五、参考文献

| 编号 | 文献 |
|---|---|
| Sabins 1999 | Sabins, F. F. (1999). Remote sensing for mineral exploration. *Ore Geology Reviews*, 14(3-4), 157-183. |
| Crósta & Moore 1989 | Crósta, A. P., & Moore, J. McM. (1989). Enhancement of Landsat Thematic Mapper imagery for residual soil mapping in SW Minas Gerais State, Brazil. *Proc. 7th Thematic Conf. Remote Sensing for Exploration Geology*. |
| Loughlin 1991 | Loughlin, W. P. (1991). Principal component analysis for alteration mapping. *Photogrammetric Engineering & Remote Sensing*, 57(9), 1163-1169. |
| Rowan & Mars 2003 | Rowan, L. C., & Mars, J. C. (2003). Lithologic mapping in the Mountain Pass, California area using ASTER data. *Remote Sensing of Environment*, 84(3), 350-366. |
| Mars & Rowan 2006 | Mars, J. C., & Rowan, L. C. (2006). Regional mapping of phyllic- and argillic-altered rocks in the Zagros magmatic arc, Iran, using ASTER data and logical operator algorithms. *Geosphere*, 2(3), 161-186. |
| Hunt 1977 | Hunt, G. R. (1977). Spectral signatures of particulate minerals in the visible and near infrared. *Geophysics*, 42(3), 501-513. |
| Clark et al. (USGS Spectral Library) | Clark, R. N., et al. USGS Digital Spectral Library splib06/07. https://www.usgs.gov/labs/spectroscopy-lab/science/spectral-library |
| Saunders et al. 1999 | Saunders, D. F., Burson, K. R., & Thompson, C. K. (1999). Model for hydrocarbon microseepage and related near-surface alteration. *AAPG Bulletin*, 83(1), 170-185. |
| van der Meer et al. 2012 | van der Meer, F. D., et al. (2012). Multi- and hyperspectral geologic remote sensing: A review. *International Journal of Applied Earth Observation and Geoinformation*, 14(1), 112-128. |
| Crósta et al. 2003 | Crósta, A. P., De Souza Filho, C. R., Azevedo, F., & Brodie, C. (2003). Targeting key alteration minerals in epithermal deposits in Patagonia, Argentina, using ASTER imagery and PCA. *International Journal of Remote Sensing*, 24(21), 4233-4240. |
| Cudahy 1997 | Cudahy, T. J., et al. (1997). Mapping porphyry-skarn alteration at Yerington, Nevada, using airborne hyperspectral VNIR-SWIR data. *Aust. CSIRO Report*. |
| Mars 2018 | Mars, J. C. (2018). Mineral and lithologic mapping capability of WorldView 3 data at Mountain Pass, California, using true- and false-color composite images, band ratios, and logical operator algorithms. *Economic Geology*, 113(7), 1587-1601. |
| Pour & Hashim 2014 | Pour, A. B., & Hashim, M. (2014). ASTER, ALI and Hyperion sensors data for lithological mapping and ore minerals exploration. *Springer Plus*, 3, 130. |
| Pour & Hashim 2015 | Pour, A. B., & Hashim, M. (2015). Hydrothermal alteration mapping from Landsat-8 data, Sar Cheshmeh copper mining district, south-eastern Islamic Republic of Iran. *Journal of Taibah University for Science*, 9(2), 155-166. |
| Yousefi & Carranza 2015 | Yousefi, M., & Carranza, E. J. M. (2015). Geometric average of spatial evidence data layers: A GIS-based multi-criteria decision-making approach to mineral prospectivity mapping. *Computers & Geosciences*, 83, 72-79. |
| Rockwell 2013 | Rockwell, B. W. (2013). Automated mapping of mineral groups and green vegetation from Landsat Thematic Mapper imagery with an example from the San Juan Mountains, Colorado. *USGS Scientific Investigations Map 3252*. |
| Bedini 2009 | Bedini, E. (2009). Mapping lithology of the Sarfartoq carbonatite complex, southern West Greenland, using HyMap imaging spectrometer data. *Remote Sensing of Environment*, 113(6), 1208-1219. |
| Noomen et al. 2008 | Noomen, M. F., et al. (2008). Hyperspectral indices for detecting changes in canopy reflectance as a result of underground natural gas leakage. *International Journal of Remote Sensing*, 29(20), 5987-6008. |
| Cloutis 1989 | Cloutis, E. A. (1989). Spectral reflectance properties of hydrocarbons: remote-sensing implications. *Science*, 245(4914), 165-168. |
| Yang et al. 2000 | Yang, H., et al. (2000). Geochemical and mineralogical anomalies in soils above hydrocarbon-bearing sediments. *AAPG Bulletin*. |
| Khan & Jacobson 2008 | Khan, S. D., & Jacobson, S. (2008). Remote sensing in unconventional resource exploration. *AAPG Memoir*. |
| Pirajno 2009 | Pirajno, F. (2009). *Hydrothermal Processes and Mineral Systems*. Springer. |
| Skirrow 2010 | Skirrow, R. G. (2010). "Hematite-group" IOCG ± U ore systems: Tectonic settings, hydrothermal characteristics, and Cu-Au and U mineralizing processes. *Geol. Assoc. Canada Short Course Notes*, 20. |
| Leach et al. 2010 | Leach, D. L., Bradley, D. C., Huston, D., Pisarevsky, S. A., Taylor, R. D., & Gardoll, S. J. (2010). Sediment-hosted lead-zinc deposits in Earth history. *Economic Geology*, 105(3), 593-625. |
| Slack 1996 | Slack, J. F. (1996). Tourmaline associations with hydrothermal ore deposits. *Reviews in Mineralogy*, 33, 559-643. |
| Hunt & Salisbury 1976 | Hunt, G. R., & Salisbury, J. W. (1976). Visible and near-infrared spectra of minerals and rocks: XI. Sedimentary rocks. *Modern Geology*, 5, 211-217. |
| Cui et al. 2019 | Cui, J., Yan, B., Dong, X., et al. (2019). Mineralogy and geochemistry of Mn-rich skarns in the Shizhuyuan polymetallic deposit, South China. *Ore Geology Reviews*. |
| Hofstra & Cline 2000 | Hofstra, A. H., & Cline, J. S. (2000). Characteristics and models for Carlin-type gold deposits. *Reviews in Economic Geology*, 13, 163-220. |
| Rowan & Bowers 1995 | Rowan, L. C., & Bowers, T. L. (1995). Analysis of linear features mapped in Landsat thematic mapper and side-looking airborne radar images of the Reno 1° x 2° quadrangle, Nevada and California. *USGS Professional Paper 1538-G*. |
| Sillitoe 2010 | Sillitoe, R. H. (2010). Porphyry copper systems. *Economic Geology*, 105(1), 3-41. |
| Hauff 2008 | Hauff, P. L. (2008). An overview of VNIR-SWIR field spectroscopy as applied to mineral exploration. *Spectral International Inc. Technical Report*. |
| Hewett 1968 | Hewett, D. F. (1968). Silver in veins of hypogene manganese oxides. *USGS Circular*, 553. |
| Rajendran & Nasir 2014 | Rajendran, S., & Nasir, S. (2014). Hydrothermal altered serpentinized zone and a study of Ni-magnesioferrite-magnetite-awaruite occurrences in Wadi Hibi, Northern Oman Mountain. *Ore Geology Reviews*, 62, 211-226. |
| Rajendran et al. 2011 | Rajendran, S., Al-Khirbash, S., Pracejus, B., Nasir, S., Al-Abri, A. H., Kusky, T. M., & Ghulam, A. (2011). ASTER detection of chromite bearing mineralized zones in Semail Ophiolite Massifs of the northern Oman Mountains. *Ore Geology Reviews*, 44, 121-135. |
| Carten et al. 1993 | Carten, R. B., White, W. H., & Stein, H. J. (1993). High-grade granite-related molybdenum systems: Classification and origin. *Geological Association of Canada Special Paper*, 40, 521-554. |
| Meinert et al. 2005 | Meinert, L. D., Dipple, G. M., & Nicolescu, S. (2005). World skarn deposits. *Economic Geology 100th Anniversary Volume*, 299-336. |
| Černý et al. 2005 | Černý, P., Blevin, P. L., Cuney, M., & London, D. (2005). Granite-related ore deposits. *Economic Geology 100th Anniversary Volume*, 337-370. |
| Hitzman et al. 2005 | Hitzman, M. W., Selley, D., & Bull, S. (2005). Formation of sedimentary rock-hosted stratiform copper deposits through Earth history. *Economic Geology*, 100(4), 609-642. |
| Berger 2011 | Berger, V. I., Singer, D. A., Bliss, J. D., & Moring, B. C. (2011). Ni-Co laterite deposits of the world—database and grade and tonnage models. *USGS Open-File Report 2011-1058*. |
| Verplanck et al. 2014 | Verplanck, P. L., Mariano, A. N., & Mariano, A. Jr. (2014). Rare earth element ore geology of carbonatites. *Reviews in Economic Geology*, 18, 5-32. |
| Mars & Rowan 2011 | Mars, J. C., & Rowan, L. C. (2011). ASTER spectral analysis and lithologic mapping of the Khanneshin carbonatite volcano, Afghanistan. *Geosphere*, 7(1), 276-289. |
| Bao & Zhao 2008 | Bao, Z., & Zhao, Z. (2008). Geochemistry of mineralization with exchangeable REY in the weathering crusts of granitic rocks in South China. *Ore Geology Reviews*, 33(3-4), 519-535. |
| Cardoso-Fernandes et al. 2019 | Cardoso-Fernandes, J., Teodoro, A. C., & Lima, A. (2019). Remote sensing data in lithium (Li) exploration: A new approach for the detection of Li-bearing pegmatites. *International Journal of Applied Earth Observation and Geoinformation*, 76, 10-25. |
| Cardoso-Fernandes et al. 2021 | Cardoso-Fernandes, J., et al. (2021). Tools for remote exploration: a lithium (Li) dedicated spectral library of the Fregeneda-Almendra aplite-pegmatite field. *Data*, 6(3), 33. |
| Kesler et al. 2012 | Kesler, S. E., Gruber, P. W., Medina, P. A., Keoleian, G. A., Everson, M. P., & Wallington, T. J. (2012). Global lithium resources: Relative importance of pegmatite, brine and other deposits. *Ore Geology Reviews*, 48, 55-69. |
| Jefferson et al. 2007 | Jefferson, C. W., et al. (2007). Unconformity-associated uranium deposits of the Athabasca Basin, Saskatchewan and Alberta. *Geological Association of Canada Mineral Deposits Division Special Publication*, 5, 273-305. |
| Cuney et al. 2012 | Cuney, M., Emetz, A., Mercadier, J., Mykchaylov, V., Shunko, V., & Yuslenko, A. (2012). Uranium deposits associated with Na-metasomatism from central Ukraine: A review of some of the major deposits and genetic constraints. *Ore Geology Reviews*, 44, 82-106. |
| Bárdossy 1982 | Bárdossy, G. (1982). Karst Bauxites: Bauxite Deposits on Carbonate Rocks. *Developments in Economic Geology*, 14, Elsevier. |
| Mitchell 1986 | Mitchell, R. H. (1986). *Kimberlites: Mineralogy, Geochemistry, and Petrology*. Springer. |

---

## 六、使用建议

1. **优先用 ASTER**:6 个 SWIR 波段是黏土族 / 碳酸盐 / Mg-OH 细分的"金标准"。Landsat 8 / Sentinel-2 在 SWIR 只有 1–2 个宽带,无法细分(只能粗判"有 vs 无 Al-OH 吸收")。
2. **波段比值法 vs Crósta PCA 的选择**:
   - **比值法** 公式固定、可解释强、对小区域稳定;缺点是受地形 / 大气残留 / 暗化背景影响大,阈值难定。
   - **Crósta PCA** 在 ROI 内自适应,能压制岩石本底和地形阴影,信噪比高;缺点是单 ROI 像素数过少(< 几千)时 PC 不稳定。
   - **生产实践**:两种都跑,取交集 → 高置信靶区(本系统综合判断 tab 即此用法)。
3. **新增矿种 / 蚀变入库时**:
   - 元素 → 矿物 → 中心吸收波长 → 卫星波段命中 → 比值公式 / PCA 子集
   - 至少要有 2 篇主流文献支持(如本表列出的文献)
   - 入库后在已知矿区跑 1–2 次回归测试再上线
4. **Sentinel-2 别期望太多**:它的 SWIR 只能做"是否有黏土"二分,做不了"绢云母 vs 高岭石"。如果项目仅有 Sentinel-2 数据,矿床类型识别只能停留在很粗的"含羟基蚀变带"层面。
5. **油气微渗漏**:近地表蚀变弱、面积大,**单指数易误判**;红层褪色 + 次生碳酸盐化 + 复合烃指数 + 植被红边胁迫 NDRE 四件套要一起跑,加已知井位标定后再圈靶。
