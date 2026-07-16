# -*- coding: utf-8 -*-
"""V1.0 冻结验证一键脚本:
① 默认参数张强基准 (冒烟基准重标用)  ② 冻结参数张强终点+曲线特征 (论文素材04回填)
③ 冻结参数固化 data/calibrated/v1_frozen_params.json
结果存 data/results/freeze_v1_output.txt"""
import sys, io, os, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)   # scripts/ 上一级 = 项目根
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
from model.engine import Engine

lines = []
def emit(s=""):
    print(s); lines.append(s)

FINAL = {'alpha_steel': 0.77182, 'alpha_slag': 1.16068, 'G_CO0': 1.22544, 'gamma': 0.44005,
         'eta_post': 0.10035, 'Q_rad_W': 2224473.16, 'A_emulsion': 317.5931, 'C_direct': 0.24952,
         'eta_O2_late': 0.45, 'scrap_tau_heat': 209.06472, 'scrap_heat_scale': 0.60007,
         'aO_cap_CO_mult': 2.04938}
z = {'铁水质量':170.0,'铁水温度':1276.0,'w([C])':4.48,'w([Si])':0.34,'w([Mn])':0.19,'w([P])':0.117,
 'w([S])':0.01,'w([V])':0.0,'w([Ti])':0.0,'废钢加入量':17.0,'供氧流量':850.0,'底吹氩气流量':11.5,
 '炉渣碱度':3.0,'MgO目标含量':8.0,'终点碳目标':0.04,'终点温度目标':1650.0,
 '枪位高度':1.5,'喷嘴直径':0.05,'喷嘴夹角':12.0,'熔池深度':1.8}

r0 = Engine().run(z)
emit("[基准] 默认参数张强炉: C=%.4f T=%.1f  (冒烟[1]重标为此值)" % (r0['终点C'], r0['终点T']))

r = Engine(params=FINAL).run(z); ts = r['时序数据']
Tc = ts['T_C']; tmin = int(np.argmin(Tc[:120]))
feo = np.array(ts['w_FeO']); pk = int(np.argmax(feo[:len(feo)//2]))
vl = int(np.argmin(feo[pk:len(feo)*4//5])) + pk
emit("[V1.0] 冻结参数张强炉: C=%.4f T=%.1f 时长%.0fs(%.1fmin)" % (
    r['终点C'], r['终点T'], ts['time_s'][-1], ts['time_s'][-1]/60))
emit("  终渣: FeO=%.1f CaO=%.1f SiO2=%.1f MgO=%.1f MnO=%.2f P2O5=%.2f" % (
    ts['w_FeO'][-1], ts['w_CaO'][-1], ts['w_SiO2'][-1], ts['w_MgO'][-1], ts['w_MnO'][-1], ts['w_P2O5'][-1]))
emit("  终点钢: Mn=%.4f P=%.4f" % (ts['w_Mn'][-1], ts['w_P'][-1]))
emit("  曲线: 微降 起%.1f 最低%.1f@%.0fs | FeO 峰%.1f@%.0fs 谷%.1f@%.0fs | Si归零 %.0fs" % (
    Tc[0], Tc[tmin], ts['time_s'][tmin], feo[pk], ts['time_s'][pk], feo[vl], ts['time_s'][vl],
    next((ts['time_s'][i] for i in range(len(ts['w_Si'])) if ts['w_Si'][i] < 0.01), -1)))
emit("  氧账本: %s" % r['氧账本'])

with open(os.path.join(ROOT, 'data', 'calibrated', 'v1_frozen_params.json'), 'w', encoding='utf-8') as f:
    json.dump(FINAL, f, ensure_ascii=False, indent=2)
emit("[固化] 冻结参数已写入 data/calibrated/v1_frozen_params.json")

out = os.path.join(ROOT, 'data', 'results', 'freeze_v1_output.txt')
with open(out, 'w', encoding='utf-8') as f:
    f.write("\n".join(lines) + "\n")
print("\n[结果已保存到 %s]" % out)
