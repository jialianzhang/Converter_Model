# 基于耦合反应机理的转炉冶炼模型

在静态模型基础上,基于耦合反应机理建立的转炉(BOF)冶炼过程模型:静态模型按物料平衡+热量平衡反算辅料/耗氧/合理废钢量,动态模型以单一界面氧活度为纽带、按 13 步时间循环推进,输出熔池温度、钢水成分、炉渣成分随吹炼时间的全程变化与终点预测。架构对标张强(2024,《钢铁》)两段式 + Kitamura(2009, ISIJ)传质/溶解动力学。

**当前版本:V1.0(已冻结,git tag `v1.0`)。** `Engine()` 不带参数即是生产模型(12 个设备级标定参数已烘入 `config.CALIBRATION_DEFAULTS`)。任何模型改进走新版本,不回改 v1.0。

## 快速开始

```bash
python main.py                     # 单炉演示(张强炉工况)
python main.py data/raw/xxx.csv    # 批量预测(CSV 格式见 data/raw/heats_template.csv)
python -m pytest tests/ -q         # 53 个单元测试
python scripts/_calib_smoke.py     # 回归守卫(基准: 张强炉 C=0.0292 / T=1700.3)
```

所有脚本从**项目根目录**运行。核心接口:

```python
from model.engine import Engine
r = Engine().run(inputs)   # → {终点C, 终点T, 时序数据(含FeO三通道/a_O_star/全成分), 氧账本, 热账本_MJ}
# inputs['终止模式']='目标碳' → 使用口径(吹到目标碳停, 张强同口径)
# 默认(提供 实际耗氧量Nm3 时) → 氧锚定口径(吹完实际氧量读终点, 检验用)

from model.static_model import StaticModel
s = StaticModel().calculate(inputs)  # → 配料/耗氧/合理废钢/热平衡明细/计算vs实际校验
```

输入字典键见 `iomodules/input_reader.py`。

## 目录地图

| 目录/文件 | 内容 |
|---|---|
| `config.py` | 全局常数 + **V1.0 冻结标定参数**(CALIBRATION_DEFAULTS)与物理范围 |
| `main.py` | 入口:单炉 / 批量 / 标定 |
| `model/` | 动态引擎(engine.py 为核心 13 步循环)、静态模型、传质/混合能/溶解/热平衡子模块 |
| `thermo/` | 热力学数据与模型:平衡常数(JSPS)、Wagner 活度、渣相(Riboud 黏度/γ_FeO/L_P)、离子理论 |
| `calibration/` | 参数标定:目标函数(objective)与优化器(optimizer, Nelder-Mead 有界 + 70/30 划分) |
| `iomodules/` | CSV 读入(实测列解析)与结果输出 |
| `tests/` | 53 个单元测试(pytest) |
| `scripts/` | 工作流脚本,见下表 |
| `data/` | `raw/` 生产数据(**专有,勿入库**,.gitignore 已挡);`calibrated/` 标定产物;`results/` 运行输出与诊断存档 |
| `history/` | **论文素材 00~05**(总览/静态/动态公式/标定泛化/结果/讨论局限)+ `配图/`(10 张 300dpi 图 + 图注正文文字) |
| `modification/` | **修改编年史**`转炉模型V4修改方案.md`(13.x 节,V1→冻结全过程)+ 早期审阅建议 |
| `formula/` | 公式文档(注意:均早于 V4~V1.0 修正,**以代码与素材 02 为准**;旧版在 `archive/`) |
| `docs/archive/` | V2 时代设计/评审文档(历史存档) |
| `papers/` | 参考文献(张强 2024、Kitamura 2009 等) |

## scripts/ 工作流脚本

| 脚本 | 用途 | 时机 |
|---|---|---|
| `_calib_smoke.py` | 回归守卫:冻结基准复现 + 参数敏感性方向 + 标定链路冒烟 | 任何改动后必跑 |
| `_paper_stats.py` | **论文数字统一快照**(泛化表/渣组分/代表炉曲线特征/两本账)→ data/results/paper_stats_v1.txt | 论文数字被质疑时一键复现 |
| `_gen_figures.py` | 10 张论文配图重生成 → history/配图/ | 图面调整后 |
| `_run_final_calibration.py` | 最终标定(目标碳框架,eta_O2_late 固定 0.45) | 换炉役/喷头/原料结构大变时重标 |
| `_run_calibration.py` | 首轮标定脚本(氧锚定框架);`check` 模式=纯数据体检 | 新数据批次先跑 check |
| `_plot_aO.py` | 界面氧活度三面板诊断图 | 按需 |
| `_freeze_v1.py` | V1.0 冻结一键验证+参数固化(历史一次性) | 存档 |

## 三个"以哪为准"(防文档漂移)

1. **公式**:以代码为准(用户规则:先改代码后改文档);成文版以 `history/论文素材_02` 为准。`formula/` 各版均为 V4 修正前的旧文档。
2. **数字**:以 `scripts/_paper_stats.py` 冻结实算快照为准(`data/results/paper_stats_v1.txt`)。标定报告用全精度参数评估,与冻结舍入值有 ≤0.4°C 漂移,勿混用。
3. **历史**:以 `modification/转炉模型V4修改方案.md` 编年史为准(每次改动的动机、验证与负结果都在)。

## 论文口径红线

- 终点**碳**命中率由目标碳终止条件构造(张强同口径),不作自由预测证据;**温度/锰/磷**才是自由预测。
- 模型 a_O* 为**界面**氧活度,副枪 TSO 实测为**本体**溶解氧,存在物理必然的系统差(约 3.1×)。
- 所有曲线与数字均为公式实算,禁止为凑参考曲线调参。

## 数据安全

`data/raw/heats_100.csv`、`data/raw/heats.csv`、`data/rael/`(全炉役台账)为**企业专有生产数据**,已列入 .gitignore,任何情况下不得提交入库;仓库内仅保留脱敏模板 `data/raw/heats_template.csv`。
