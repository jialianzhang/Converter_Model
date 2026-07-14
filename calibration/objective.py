"""参数标定目标函数

核心原则: 一组共享参数联合最小化 N 炉的预测误差 —— 标定的是设备/工艺级常数
(二次燃烧率、散热、乳化面积等), 不是逐炉拟合; 炉次间差异由每炉输入变量承载。

多指标联合约束: C/T 必标, FeO 可选 —— V4.2 已证明终点 FeO 与 eta_post 存在
补偿关系(FeO 过冲放热回填热缺口), 只标 C/T 会让两个误差互相掩护。
"""
import numpy as np
from model.engine import Engine

# 参数向量顺序 (与 config.CALIBRATION_DEFAULTS / CALIBRATION_BOUNDS 的键对应)
PARAM_NAMES = ['alpha_steel', 'alpha_slag', 'G_CO0', 'gamma', 'eta_post',
               'Q_rad_W', 'A_emulsion', 'C_direct', 'eta_O2_late', 'scrap_tau_heat',
               'scrap_heat_scale', 'aO_cap_CO_mult']

# 误差归一尺度 (使不同量纲的指标在目标函数中权重可比):
# ΔC 按 0.01% (命中带宽的一半), ΔT 按 10°C (命中带宽), ΔFeO 按 5%,
# Δa_O 按 0.02 (≈200ppm 溶解氧), Δ出钢量 按 2 t, ΔMn 按 0.02%, ΔP 按 0.005%
SCALE_C, SCALE_T, SCALE_FEO, SCALE_AO, SCALE_W = 0.01, 10.0, 5.0, 0.02, 2.0
SCALE_MN, SCALE_P = 0.02, 0.005


def vector_to_params(x) -> dict:
    return dict(zip(PARAM_NAMES, [float(v) for v in x]))


def params_to_vector(p: dict):
    return np.array([p[k] for k in PARAM_NAMES], dtype=float)


def _eval_one(args):
    """单炉评估 (可作 multiprocessing worker); 单炉仿真失败返回 None, 由上层计罚

    framing='oxygen': 氧锚定 (吹真实氧量, 碳温皆自由 —— 最严苛检验)
    framing='target_c': 目标碳终止 (张强口径/使用框架; 碳按实测碳收口=构造性命中,
                        温度/Mn/P 为真预测; 实际耗氧量仍锚定加料窗口)"""
    params, heat, framing = args
    try:
        h = heat
        if framing == 'target_c':
            h = dict(heat)
            h['终止模式'] = '目标碳'
            h['终点碳目标'] = heat['终点C_actual']
        engine = Engine(params=params)
        r = engine.run(h)
        ts = r['时序数据']
        dC = r['终点C'] - heat['终点C_actual']
        dT = r['终点T'] - heat['终点T_actual']
        dFeO = None
        if heat.get('终点FeO_actual') is not None:
            dFeO = ts['w_FeO'][-1] - heat['终点FeO_actual']
        daO = None
        if heat.get('终点O_actual_ppm') is not None:
            daO = ts['a_O_star'][-1] - heat['终点O_actual_ppm'] / 1e4
        dW = None
        if heat.get('出钢量_actual') is not None:
            dW = ts['m_steel_kg'][-1] / 1000.0 - heat['出钢量_actual']
        dMn = None
        if heat.get('终点Mn_actual') is not None:
            dMn = ts['w_Mn'][-1] - heat['终点Mn_actual']
        dP = None
        if heat.get('终点P_actual') is not None:
            dP = ts['w_P'][-1] - heat['终点P_actual']
        return (dC, dT, dFeO, daO, dW, dMn, dP)
    except Exception:
        return None


def evaluate_heats(params: dict, heats: list, pool=None, framing='oxygen') -> list:
    """用同一组参数跑全部炉次, 返回逐炉误差 [(dC,dT,dFeO,daO,dW,dMn,dP)|None, ...]
    可选指标对应列缺失时为 None。pool: multiprocessing.Pool (100炉×数百次评估时必要)"""
    tasks = [(params, h, framing) for h in heats]
    if pool is not None:
        return pool.map(_eval_one, tasks)
    return [_eval_one(t) for t in tasks]


def calibration_objective(x, heats, w_C=1.0, w_T=1.0, w_FeO=1.0,
                          w_aO=0.5, w_W=0.5, w_Mn=0.0, w_P=0.0,
                          pool=None, framing='oxygen') -> float:
    """加权归一均方误差 (无量纲); 可选指标仅在对应实测列存在时参与"""
    rows = evaluate_heats(vector_to_params(x), heats, pool=pool, framing=framing)
    loss = 0.0
    for row in rows:
        if row is None:
            loss += 1e4   # 单炉仿真失败罚项
            continue
        dC, dT, dFeO, daO, dW, dMn, dP = row
        loss += w_C * (dC / SCALE_C) ** 2 + w_T * (dT / SCALE_T) ** 2
        if dFeO is not None:
            loss += w_FeO * (dFeO / SCALE_FEO) ** 2
        if daO is not None:
            loss += w_aO * (daO / SCALE_AO) ** 2
        if dW is not None:
            loss += w_W * (dW / SCALE_W) ** 2
        if dMn is not None:
            loss += w_Mn * (dMn / SCALE_MN) ** 2
        if dP is not None:
            loss += w_P * (dP / SCALE_P) ** 2
    return loss / max(len(rows), 1)


def hit_rates(rows) -> dict:
    """按张强论文口径统计命中率: |ΔC|≤0.02%, |ΔT|≤10°C (自动跳过失败炉次)"""
    rows = [r for r in rows if r is not None]
    n = max(len(rows), 1)
    hC = sum(1 for row in rows if abs(row[0]) <= 0.02) / n
    hT = sum(1 for row in rows if abs(row[1]) <= 10.0) / n
    out = {'命中率_C±0.02': round(hC, 3), '命中率_T±10': round(hT, 3),
           'MAE_C': round(float(np.mean([abs(r[0]) for r in rows])), 4),
           'MAE_T': round(float(np.mean([abs(r[1]) for r in rows])), 1)}
    feo_rows = [r[2] for r in rows if r[2] is not None]
    if feo_rows:
        out['MAE_FeO'] = round(float(np.mean([abs(v) for v in feo_rows])), 1)
    ao_rows = [r[3] for r in rows if r[3] is not None]
    if ao_rows:
        out['MAE_aO'] = round(float(np.mean([abs(v) for v in ao_rows])), 4)
    w_rows = [r[4] for r in rows if r[4] is not None]
    if w_rows:
        out['MAE_出钢量t'] = round(float(np.mean([abs(v) for v in w_rows])), 2)
    mn_rows = [r[5] for r in rows if len(r) > 5 and r[5] is not None]
    if mn_rows:
        out['MAE_Mn'] = round(float(np.mean([abs(v) for v in mn_rows])), 4)
        out['命中率_Mn±0.02'] = round(sum(1 for v in mn_rows if abs(v) <= 0.02) / len(mn_rows), 3)
    p_rows = [r[6] for r in rows if len(r) > 6 and r[6] is not None]
    if p_rows:
        out['MAE_P'] = round(float(np.mean([abs(v) for v in p_rows])), 4)
        out['命中率_P±0.005'] = round(sum(1 for v in p_rows if abs(v) <= 0.005) / len(p_rows), 3)
    return out
