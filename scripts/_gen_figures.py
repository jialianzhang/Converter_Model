# -*- coding: utf-8 -*-
"""V1.0 论文配图批量生成 → history/配图/fig01~10 (300dpi PNG)
风格: 中文标签(微软雅黑), Okabe-Ito 色盲安全配色, 统一字号, 期刊单/双栏尺寸。
数据: V1.0 冻结模型 (Engine() 默认即冻结参数), 100 炉目标碳框架。"""
import sys, io, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)   # scripts/ 上一级 = 项目根
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DengXian', 'Arial'],
    'axes.unicode_minus': False,
    'font.size': 9, 'axes.labelsize': 9.5, 'axes.titlesize': 10.5,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'legend.fontsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.8, 'lines.linewidth': 1.6,
    'figure.dpi': 110, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})
import matplotlib.pyplot as plt
from model.engine import Engine
from thermo.wagner import calc_lgf_i
from iomodules.input_reader import InputReader

# Okabe-Ito 色盲安全
C_BLUE, C_ORANGE, C_GREEN = '#0072B2', '#E69F00', '#009E73'
C_RED, C_PURPLE, C_SKY, C_GREY = '#D55E00', '#CC79A7', '#56B4E9', '#8A8A8A'

OUT = os.path.join(ROOT, 'history', '配图')
os.makedirs(OUT, exist_ok=True)
def save(fig, name):
    p = os.path.join(OUT, name); fig.savefig(p); plt.close(fig); print('  已生成', name)

# ---------- 数据准备 ----------
heats = [InputReader.parse_single_heat(x) for x in InputReader.read_csv(os.path.join(ROOT, 'data', 'raw', 'heats_100.csv'))]
for h in heats:
    if h.get('实际烧结矿kg') is None: h['实际烧结矿kg'] = 0.0
# 代表炉次: 各关键量最接近中位数
keys = ['铁水质量', '废钢加入量', '供氧流量', '终点C_actual', '终点T_actual']
med = {k: np.median([h[k] for h in heats]) for k in keys}
rep = min(heats, key=lambda h: sum(((h[k]-med[k])/max(abs(med[k]),1e-6))**2 for k in keys))
eng = Engine()   # V1.0 冻结模型
def run_tc(h):
    h2 = dict(h); h2['终止模式'] = '目标碳'; h2['终点碳目标'] = h['终点C_actual']
    return eng.run(h2)
print('代表炉: 铁水%.1ft 废钢%.1ft 实测C%.3f T%.0f' % (rep['铁水质量'], rep['废钢加入量'], rep['终点C_actual'], rep['终点T_actual']))
r_rep = run_tc(rep); ts = r_rep['时序数据']
t_min = np.array(ts['time_s'])/60.0

print('批量跑 100 炉 (散点/直方图用)...')
pred = []
for h in heats:
    r = run_tc(h); s = r['时序数据']
    pred.append(dict(T=r['终点T'], Mn=s['w_Mn'][-1], P=s['w_P'][-1], FeO=s['w_FeO'][-1],
                     CaO=s['w_CaO'][-1], SiO2=s['w_SiO2'][-1], MgO=s['w_MgO'][-1],
                     MnO=s['w_MnO'][-1], P2O5=s['w_P2O5'][-1]))
actT = np.array([h['终点T_actual'] for h in heats]); modT = np.array([p['T'] for p in pred])
actMn = np.array([h['终点Mn_actual'] for h in heats]); modMn = np.array([p['Mn'] for p in pred])
# 训练/验证划分 (与 Calibrator 同种子)
idx = np.random.default_rng(42).permutation(len(heats)); tr, va = idx[:70], idx[70:]

# ---------- 图1: 单炉全程曲线三联 (对照张强图4) ----------
fig, ax = plt.subplots(1, 3, figsize=(10.8, 3.1))
ax[0].plot(t_min, ts['T_C'], color=C_RED)
ax[0].set_xlabel('吹炼时间/min'); ax[0].set_ylabel('熔池温度/℃'); ax[0].set_title('(a) 熔池温度')
ax[0].axhline(rep['终点T_actual'], ls='--', lw=0.9, color=C_GREY)
ax[0].annotate('实测终点 %d ℃' % rep['终点T_actual'], xy=(0.98, rep['终点T_actual']),
               xycoords=('axes fraction', 'data'), ha='right', va='bottom', fontsize=8, color=C_GREY)
a1 = ax[1]; a1b = a1.twinx(); a1b.spines['right'].set_visible(True)
a1.plot(t_min, ts['w_C'], color=C_BLUE, label='C')
a1b.plot(t_min, ts['w_Si'], color=C_ORANGE, label='Si')
a1b.plot(t_min, ts['w_Mn'], color=C_GREEN, label='Mn')
a1b.plot(t_min, ts['w_P'], color=C_PURPLE, label='P')
a1.set_xlabel('吹炼时间/min'); a1.set_ylabel('w(C)/%'); a1b.set_ylabel('w(Si, Mn, P)/%')
a1.set_title('(b) 熔池成分'); h1,l1 = a1.get_legend_handles_labels(); h2,l2 = a1b.get_legend_handles_labels()
a1.legend(h1+h2, l1+l2, ncol=2, frameon=False, loc='upper right')
ax[2].plot(t_min, ts['w_FeO'], color=C_RED, label='FeO')
ax[2].plot(t_min, ts['w_CaO'], color=C_BLUE, label='CaO')
ax[2].plot(t_min, ts['w_SiO2'], color=C_ORANGE, label='SiO$_2$')
ax[2].plot(t_min, ts['w_MgO'], color=C_GREEN, label='MgO')
ax[2].plot(t_min, ts['w_MnO'], color=C_PURPLE, label='MnO')
ax[2].set_xlabel('吹炼时间/min'); ax[2].set_ylabel('炉渣组分/%'); ax[2].set_title('(c) 炉渣成分')
ax[2].legend(ncol=2, frameon=False)
fig.tight_layout(); save(fig, 'fig01_单炉全程曲线三联.png')

# ---------- 图2: FeO 三通道源汇分解 ----------
tarr = np.array(ts['time_s']); tend = tarr[-1]
gen = np.array(ts['feo_gen_impact']); ifc = np.array(ts['feo_net_iface']); dep = np.array(ts['feo_cons_deP'])
labels = ['前期\n(0~33%)', '中期\n(33~66%)', '末期\n(66~100%)']
G, I, D = [], [], []
for a, b in [(0, .33), (.33, .66), (.66, 1.001)]:
    m = (tarr >= a*tend) & (tarr < b*tend)
    G.append(gen[m].sum()/1000); I.append(ifc[m].sum()/1000); D.append(dep[m].sum()/1000)
x = np.arange(3); w = 0.27
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.2))
ax[0].bar(x-w, G, w, color=C_ORANGE, label='冲击区生成 Fe+½O$_2$→FeO')
ax[0].bar(x,   I, w, color=C_BLUE, label='界面净通量(负=被C还原)')
ax[0].bar(x+w, D, w, color=C_PURPLE, label='脱磷消耗')
net = np.array(G)+np.array(I)+np.array(D)
ax[0].plot(x, net, 'o-', color='k', lw=1.2, ms=4, label='净变化')
ax[0].axhline(0, lw=0.8, color='k'); ax[0].set_xticks(x); ax[0].set_xticklabels(labels)
ax[0].set_ylabel('FeO 质量变化/t'); ax[0].set_title('(a) FeO 三通道分期收支'); ax[0].legend(frameon=False, fontsize=7.5)
ax[1].plot(t_min, ts['w_FeO'], color=C_RED)
ax[1].set_xlabel('吹炼时间/min'); ax[1].set_ylabel('w(FeO)/%'); ax[1].set_title('(b) 渣中 FeO 全程曲线')
for frac, lab in [(0.33, ''), (0.66, '')]:
    ax[1].axvline(frac*tend/60, ls=':', lw=0.9, color=C_GREY)
fig.tight_layout(); save(fig, 'fig02_FeO三通道源汇分解.png')

# ---------- 图3: 界面氧活度全程 ----------
TK = np.array(ts['T_K']); aO = np.array(ts['a_O_star'])
ppm, coeq = [], []
for i in range(len(tarr)):
    comp = {'C': ts['w_C'][i], 'Si': ts['w_Si'][i], 'Mn': ts['w_Mn'][i], 'P': ts['w_P'][i], 'S': ts['w_S'][i]}
    fO = 10.0**calc_lgf_i('O', comp, TK[i])
    ppm.append(aO[i]/max(fO, 1e-6)*1e4)
    coeq.append(0.0025/max(ts['w_C'][i], 1e-4)*1e4)
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.2))
a0 = ax[0]; a0.plot(t_min, aO, color=C_BLUE); a0.set_yscale('log')
a0.set_xlabel('吹炼时间/min'); a0.set_ylabel('界面氧活度 $a_\\mathrm{O}^*$'); a0.set_title('(a) 界面氧活度(对数)')
a0b = a0.twinx(); a0b.spines['right'].set_visible(True)
a0b.plot(t_min, ts['w_C'], color=C_GREY, ls='--', lw=1.1); a0b.set_ylabel('w(C)/%', color=C_GREY)
# 跳过开吹点火瞬态 (前 ~0.3min): 该段界面氧尖峰属首分钟瞬态(已列局限),
# 不跳则单点尖峰(~4000ppm)压扁整条曲线, 末期过冲(~2100 vs 碳氧平衡~500)读不出
ppm = np.array(ppm); coeq = np.array(coeq); mm = t_min >= 0.3
ax[1].plot(t_min[mm], ppm[mm], color=C_RED, label='模型界面氧')
ax[1].plot(t_min[mm], coeq[mm], color='k', ls='--', lw=1.2, label='碳氧平衡本体氧\n([C][O]=0.0025)')
ax[1].set_ylim(0, max(ppm[mm].max(), 1500)*1.08)
ax[1].set_xlabel('吹炼时间/min'); ax[1].set_ylabel('氧含量/ppm'); ax[1].set_title('(b) 界面氧 vs 碳氧平衡')
ax[1].legend(frameon=False)
fig.tight_layout(); save(fig, 'fig03_界面氧活度全程.png')

# ---------- 图4/5: 标定散点 T 与 Mn (训练/验证分色) ----------
def scatter_fig(act, mod, tol, unit, name, title, fmt='%g'):
    fig, ax = plt.subplots(figsize=(4.0, 3.8))
    lo = min(act.min(), mod.min()); hi = max(act.max(), mod.max())
    pad = (hi-lo)*0.08; lo, hi = lo-pad, hi+pad
    ax.fill_between([lo, hi], [lo-tol, hi-tol], [lo+tol, hi+tol], color=C_GREY, alpha=0.18,
                    label='±%s%s 命中带' % (fmt % tol, unit))
    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.9)
    ax.scatter(act[tr], mod[tr], s=22, color=C_BLUE, alpha=0.75, label='训练集(70炉)')
    ax.scatter(act[va], mod[va], s=30, color=C_ORANGE, marker='^', alpha=0.9, label='验证集(30炉)')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_xlabel('实测值'+unit); ax.set_ylabel('模型预测值'+unit); ax.set_title(title)
    ax.legend(frameon=False, loc='upper left', fontsize=7.5)
    fig.tight_layout(); save(fig, name)
scatter_fig(actT, modT, 10, '/℃', 'fig04_终点温度预测散点.png', '终点温度: 预测 vs 实测', '%d')
scatter_fig(actMn*100, modMn*100, 2, '×10$^{-2}$/%', 'fig05_终点锰预测散点.png', '终点锰: 预测 vs 实测')

# ---------- 图6: 误差分布直方图 ----------
fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.0))
dT = modT - actT; dMn = (modMn - actMn)
ax[0].hist(dT, bins=18, color=C_BLUE, alpha=0.85, edgecolor='white')
ax[0].axvline(0, color='k', lw=0.9); ax[0].axvline(-10, ls='--', lw=0.9, color=C_GREY); ax[0].axvline(10, ls='--', lw=0.9, color=C_GREY)
ax[0].set_xlabel('ΔT = 预测−实测/℃'); ax[0].set_ylabel('炉次数'); ax[0].set_title('(a) 终点温度误差分布')
ax[1].hist(dMn*100, bins=18, color=C_GREEN, alpha=0.85, edgecolor='white')
ax[1].axvline(0, color='k', lw=0.9); ax[1].axvline(-2, ls='--', lw=0.9, color=C_GREY); ax[1].axvline(2, ls='--', lw=0.9, color=C_GREY)
ax[1].set_xlabel('ΔMn ×10$^{-2}$/%'); ax[1].set_ylabel('炉次数'); ax[1].set_title('(b) 终点锰误差分布')
fig.tight_layout(); save(fig, 'fig06_误差分布直方图.png')

# ---------- 图7: 终点渣组分 vs 工业典型范围 ----------
comps = ['FeO', 'CaO', 'SiO2', 'MgO', 'MnO', 'P2O5']
disp = ['FeO', 'CaO', 'SiO$_2$', 'MgO', 'MnO', 'P$_2$O$_5$']
typ = {'FeO': (15, 30), 'CaO': (42, 50), 'SiO2': (12, 17), 'MgO': (6, 10), 'MnO': (1, 3), 'P2O5': (1, 3)}
vals = {c: np.array([p[c] for p in pred]) for c in comps}
fig, ax = plt.subplots(figsize=(6.2, 3.4))
x = np.arange(len(comps))
for i, c in enumerate(comps):
    lo, hi = typ[c]
    ax.add_patch(plt.Rectangle((i-0.32, lo), 0.64, hi-lo, color=C_GREEN, alpha=0.22, zorder=1))
    ax.errorbar(i, vals[c].mean(), yerr=[[vals[c].mean()-vals[c].min()], [vals[c].max()-vals[c].mean()]],
                fmt='o', color=C_BLUE, ms=6, capsize=4, lw=1.2, zorder=3)
ax.set_xticks(x); ax.set_xticklabels(disp)
ax.set_ylabel('终点炉渣组分/%'); ax.set_title('终点渣组分(均值±全距, 100炉) 与工业典型范围(绿带)')
fig.tight_layout(); save(fig, 'fig07_终点渣组分对照.png')

# ---------- 图8: 氧账本 ----------
ox = r_rep['氧账本']
items = [('直接氧化\n(冲击区+入铁)', ox['O_direct_mol']), ('界面注入\n(溶解氧)', ox['O_interface_inj_mol']),
         ('逸散/后燃', ox['O_escape_mol'])]
fig, ax = plt.subplots(figsize=(4.6, 3.4))
total = ox['O_blown_mol']; left = 0
colors = [C_BLUE, C_SKY, C_GREY]
for (lab, v), c in zip(items, colors):
    ax.barh(0, v/total*100, left=left, color=c, height=0.45)
    if v/total > 0.04:
        ax.text(left+v/total*50, 0, '%s\n%.1f%%' % (lab, v/total*100), ha='center', va='center', fontsize=8)
    left += v/total*100
ax.set_xlim(0, 100); ax.set_yticks([]); ax.set_xlabel('吹入氧量占比/%')
ax.set_title('氧账本: 全炉氧去向 (代表炉次, 化学利用率 %.1f%%)' % (ox['化学利用率']*100))
fig.tight_layout(); save(fig, 'fig08_氧账本.png')

# ---------- 图9: 热账本收支 ----------
hl = r_rep['热账本_MJ']
inc = [('元素氧化净热', hl['Q_oxid']+hl['Q_FeO_iface']), ('二次燃烧', hl['Q_postcomb']),
       ('C$_2$S/脱磷净热', hl['Q_C2S']+hl['Q_deP'])]
out_ = [('废钢升温熔化', hl['Q_scrap']), ('辅料分解', hl['Q_flux_decomp']),
        ('渣升温', hl['Q_slag_heat']), ('炉气显热', hl['Q_offgas']),
        ('辐射炉衬', hl['Q_rad']), ('烧结矿', hl['Q_sinter'])]
fig, ax = plt.subplots(figsize=(5.4, 3.6))
b = 0
for (lab, v), c in zip(inc, [C_RED, C_ORANGE, C_PURPLE]):
    ax.bar(0, v/1000, bottom=b, color=c, width=0.5)
    if v > sum(x[1] for x in inc)*0.05:
        ax.text(0, b+v/2000, lab+'\n%.1f GJ' % (v/1000), ha='center', va='center', fontsize=7.5)
    b += v/1000
b = 0
for (lab, v), c in zip(out_, [C_BLUE, C_SKY, C_GREEN, C_GREY, '#B0B0B0', '#D0D0D0']):
    ax.bar(1, v/1000, bottom=b, color=c, width=0.5)
    if v > sum(x[1] for x in out_)*0.05:
        ax.text(1, b+v/2000, lab+'\n%.1f GJ' % (v/1000), ha='center', va='center', fontsize=7.5)
    b += v/1000
ax.set_xticks([0, 1]); ax.set_xticklabels(['热收入', '热支出(不含钢渣升温)'])
ax.set_ylabel('热量/GJ'); ax.set_title('热账本: 全炉热收支 (代表炉次)')
fig.tight_layout(); save(fig, 'fig09_热账本收支.png')

# ---------- 图10: 泛化性对比 ----------
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.0))
mT = [np.abs(dT[tr]).mean(), np.abs(dT[va]).mean()]
mM = [np.abs(dMn[tr]).mean()*100, np.abs(dMn[va]).mean()*100]
for a, m, ylab, title, c in [(ax[0], mT, 'MAE(T)/℃', '(a) 终点温度', C_BLUE),
                             (ax[1], mM, 'MAE(Mn) ×10$^{-2}$/%', '(b) 终点锰', C_GREEN)]:
    bars = a.bar(['训练集(70)', '验证集(30)'], m, color=[c, C_ORANGE], width=0.55)
    for bar, v in zip(bars, m):
        a.text(bar.get_x()+bar.get_width()/2, v*1.02, '%.2f' % v, ha='center', fontsize=8.5)
    a.set_ylabel(ylab); a.set_title(title); a.set_ylim(0, max(m)*1.25)
fig.suptitle('泛化性: 验证集与训练集精度一致 (V1.0 冻结参数)', y=1.02, fontsize=10)
fig.tight_layout(); save(fig, 'fig10_泛化性对比.png')

print('\n全部完成 → history/配图/ (10张, 300dpi)')
