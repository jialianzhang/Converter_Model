# -*- coding: utf-8 -*-
"""论文数字统一快照 (V1.0 冻结模型, 同源性保证)

背景: 标定报告 (final_calibration_v1.txt) 的指标用优化器全精度参数计算,
而冻结进 config 的是舍入到 5 位小数的值 —— 两者在含钳位/变步长的仿真中
产生 ~0.4°C 的 MAE 漂移。论文交付的是冻结模型 Engine(), 故所有正文数字
以本脚本 (冻结参数 + 标定同款评估管线) 的输出为准。
输出: data/results/paper_stats_v1.txt
"""
import sys, io, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)   # scripts/ 上一级 = 项目根
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
from model.engine import Engine
from calibration.optimizer import Calibrator
from calibration.objective import evaluate_heats, hit_rates
from config import CALIBRATION_DEFAULTS
from thermo.wagner import calc_lgf_i
from iomodules.input_reader import InputReader

OUT = []
def P(s=''):
    print(s); OUT.append(s)

heats = [InputReader.parse_single_heat(r) for r in InputReader.read_csv(
    os.path.join(ROOT, 'data', 'raw', 'heats_100.csv'))]
for h in heats:
    if h.get('实际烧结矿kg') is None: h['实际烧结矿kg'] = 0.0

# ===== 1. 泛化表 (冻结参数, 标定同款评估: evaluate_heats + hit_rates, 同种子划分) =====
train, val = Calibrator(seed=42).split(heats)
frozen = dict(CALIBRATION_DEFAULTS)
P('=== 1. 泛化表 (V1.0 冻结参数, target_c 框架, 70/30 同种子划分) ===')
P('训练集(%d): %s' % (len(train), hit_rates(evaluate_heats(frozen, train, framing='target_c'))))
P('验证集(%d): %s' % (len(val), hit_rates(evaluate_heats(frozen, val, framing='target_c'))))

# ===== 2. 100 炉终渣组分 / 终点界面氧 / 渣量 =====
eng = Engine()
def run_tc(h):
    h2 = dict(h); h2['终止模式'] = '目标碳'; h2['终点碳目标'] = h['终点C_actual']
    return eng.run(h2)

comps = ['FeO', 'CaO', 'SiO2', 'MgO', 'MnO', 'P2O5']
acc = {c: [] for c in comps}
slag_t, iface_ppm, tso_ppm, util = [], [], [], []
for h in heats:
    r = run_tc(h); s = r['时序数据']
    for c in comps: acc[c].append(s['w_' + c][-1])
    slag_t.append(s['m_slag_kg'][-1] / 1000.0)
    comp_end = {'C': s['w_C'][-1], 'Si': s['w_Si'][-1], 'Mn': s['w_Mn'][-1],
                'P': s['w_P'][-1], 'S': s['w_S'][-1]}
    fO = 10.0 ** calc_lgf_i('O', comp_end, s['T_K'][-1])
    iface_ppm.append(s['a_O_star'][-1] / max(fO, 1e-6) * 1e4)
    if h.get('终点O_actual_ppm') is not None: tso_ppm.append(h['终点O_actual_ppm'])
    util.append(r['氧账本']['化学利用率'])
P()
P('=== 2. 100 炉终渣组分 (质量分数%, 冻结模型) ===')
for c in comps:
    a = np.array(acc[c])
    P('%-5s 均值 %5.1f  范围 [%.1f, %.1f]' % (c, a.mean(), a.min(), a.max()))
P('渣量均值 %.1f t | 终点界面氧均值 %.0f ppm (TSO 实测均值 %.0f ppm, n=%d)'
  % (np.mean(slag_t), np.mean(iface_ppm), np.mean(tso_ppm), len(tso_ppm)))
P('氧化学利用率: 均值 %.3f 范围 [%.3f, %.3f]' % (np.mean(util), np.min(util), np.max(util)))

# ===== 3. 代表炉次全套数字 (配图 fig01/02/03/08/09 的文字素材) =====
keys = ['铁水质量', '废钢加入量', '供氧流量', '终点C_actual', '终点T_actual']
med = {k: np.median([h[k] for h in heats]) for k in keys}
rep = min(heats, key=lambda h: sum(((h[k] - med[k]) / max(abs(med[k]), 1e-6)) ** 2 for k in keys))
r = run_tc(rep); s = r['时序数据']
tarr = np.array(s['time_s']); tend = tarr[-1]
P()
P('=== 3. 代表炉次 (铁水 %.1f t, 废钢 %.1f t, 废钢比 %.1f%%, 供氧 %.0f m3/min) ===' % (
    rep['铁水质量'], rep['废钢加入量'],
    rep['废钢加入量'] / (rep['铁水质量'] + rep['废钢加入量']) * 100, rep['供氧流量']))
P('终点: C %.4f (实测 %.3f) | T %.1f (实测 %.0f) | Mn %.4f (实测 %.4f) | P %.4f (实测 %.4f)' % (
    s['w_C'][-1], rep['终点C_actual'], r['终点T'], rep['终点T_actual'],
    s['w_Mn'][-1], rep.get('终点Mn_actual') or -1, s['w_P'][-1], rep.get('终点P_actual') or -1))
P('终渣: FeO %.1f CaO %.1f SiO2 %.1f MgO %.1f MnO %.2f P2O5 %.2f | 渣量 %.1f t | 时长 %.1f min' % (
    s['w_FeO'][-1], s['w_CaO'][-1], s['w_SiO2'][-1], s['w_MgO'][-1],
    s['w_MnO'][-1], s['w_P2O5'][-1], s['m_slag_kg'][-1] / 1000, tend / 60))
T = np.array(s['T_C'])
imin = int(np.argmin(T[:120]))
P('温度: 起 %.1f -> 前期最低 %.1f @%ds -> 终点 %.1f' % (T[0], T[imin], int(tarr[imin]), T[-1]))
wSi = np.array(s['w_Si']); iSi = int(np.argmax(wSi < 0.005)) if (wSi < 0.005).any() else -1
P('Si 归零: %ds (%.1f min)' % (int(tarr[iSi]), tarr[iSi] / 60))
wFeO = np.array(s['w_FeO'])
ipk = int(np.argmax(wFeO[:150])); igu = ipk + int(np.argmin(wFeO[ipk:300]))
P('FeO: 前期峰 %.1f%%@%ds -> 谷 %.1f%%@%ds -> 终点 %.1f%%' % (
    wFeO[ipk], int(tarr[ipk]), wFeO[igu], int(tarr[igu]), wFeO[-1]))
# 中期台阶(4~9 min)均值
mmid = (tarr >= 240) & (tarr <= 540)
P('FeO 中期台阶(4~9min)均值: %.1f%%' % wFeO[mmid].mean())

# FeO 三通道分期 (kg -> t)
gen = np.array(s['feo_gen_impact']); ifc = np.array(s['feo_net_iface']); dep = np.array(s['feo_cons_deP'])
P('FeO 三通道分期收支 (t):')
for (a, b), lab in zip([(0, .33), (.33, .66), (.66, 1.001)], ['前期', '中期', '末期']):
    m = (tarr >= a * tend) & (tarr < b * tend)
    P('  %s: 冲击区生成 %+.1f | 界面净 %+.1f | 脱磷 %+.2f | 净 %+.1f' % (
        lab, gen[m].sum() / 1000, ifc[m].sum() / 1000, dep[m].sum() / 1000,
        (gen[m] + ifc[m] + dep[m]).sum() / 1000))

# 界面氧 ppm 曲线特征
ppm = []
for i in range(len(tarr)):
    comp = {'C': s['w_C'][i], 'Si': s['w_Si'][i], 'Mn': s['w_Mn'][i], 'P': s['w_P'][i], 'S': s['w_S'][i]}
    fO = 10.0 ** calc_lgf_i('O', comp, s['T_K'][i])
    ppm.append(s['a_O_star'][i] / max(fO, 1e-6) * 1e4)
ppm = np.array(ppm); m03 = tarr >= 18
P('界面氧: 前中期台阶(0.3~8min)范围 %.0f~%.0f ppm | 终点 %.0f ppm | 终点碳氧平衡本体氧 %.0f ppm' % (
    ppm[(tarr >= 18) & (tarr <= 480)].min(), ppm[(tarr >= 18) & (tarr <= 480)].max(), ppm[-1],
    0.0025 / max(s['w_C'][-1], 1e-4) * 1e4))
P('氧账本: %s' % r['氧账本'])
P('热账本_MJ: %s' % r['热账本_MJ'])

# ~10 min 渣量事件定位 (石灰批次 0/3/6min 之外的质量阶跃 = 末批石灰集中溶解)
ms = np.array(s['m_slag_kg'])
dms = np.diff(ms)
big = [(int(tarr[i + 1]), dms[i]) for i in np.argsort(dms)[-6:][::-1]]
P('渣质量最大单步增量 TOP6 (时刻s, kg/步): %s' % [(t_, round(d, 0)) for t_, d in big])

with open(os.path.join(ROOT, 'data', 'results', 'paper_stats_v1.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT) + '\n')
print('\n[已保存 data/results/paper_stats_v1.txt]')
