# -*- coding: utf-8 -*-
"""标定管线冒烟测试 (非真实标定 —— 真标定需用户提供多炉实测数据)
1) 默认参数下引擎行为与 V4.2 一致 (参数化不改变结果)
2) 参数敏感性方向正确 (eta_post 增大 → 终点T升高)
3) 标定链路端到端可跑 (3 炉合成数据, 少量迭代)
结果保存到 data/results/calib_smoke_output.txt"""
import sys, io, os, copy
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # scripts/ 上一级 = 项目根
sys.path.insert(0, ROOT)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

lines = []
def emit(s=""):
    print(s)
    lines.append(s)

base = {'铁水质量':170.0,'铁水温度':1276.0,'w([C])':4.48,'w([Si])':0.34,'w([Mn])':0.19,
 'w([P])':0.117,'w([S])':0.01,'w([V])':0.0,'w([Ti])':0.0,'废钢加入量':17.0,'供氧流量':850.0,
 '底吹氩气流量':11.5,'炉渣碱度':3.0,'MgO目标含量':8.0,'冷却剂类型':'sinter','终点碳目标':0.04,
 '终点温度目标':1650.0,'枪位高度':1.5,'喷嘴直径':0.05,'喷嘴夹角':12.0,'熔池深度':1.8,
 '终点C_actual':0.0388,'终点T_actual':1655.0,'终点FeO_actual':None}

from model.engine import Engine

# 1) V1.0 冻结模型复现 (Engine()=冻结生产模型; 张强炉跨操作参考基准)
r0 = Engine().run(base)
emit("[1] V1.0冻结模型复现: C=%.4f%% T=%.1f°C (V1.0基准: 0.0292/1700.3, 张强炉跨操作)" % (r0['终点C'], r0['终点T']))
ok1 = abs(r0['终点C']-0.0292)<1e-4 and abs(r0['终点T']-1700.3)<0.5
emit("    冻结模型稳定: %s" % ("通过" if ok1 else "失败!"))

# 2) 参数敏感性方向
r_post = Engine(params={'eta_post': 0.20}).run(base)
emit("[2] eta_post 0.10→0.20: T %.1f → %.1f (应升高): %s" % (
    r0['终点T'], r_post['终点T'], "通过" if r_post['终点T'] > r0['终点T'] + 5 else "失败!"))
# 冻结默认 eta_O2_late=0.45; 升到 0.6 → 末期入铁氧保留率上升 → 终点 FeO 应升高 (单调)
r_late = Engine(params={'eta_O2_late': 0.6}).run(base)
feo0 = r0['时序数据']['w_FeO'][-1]; feo1 = r_late['时序数据']['w_FeO'][-1]
emit("[3] eta_O2_late 0.45→0.60: 终点FeO %.1f%% → %.1f%% (应升高): %s" % (
    feo0, feo1, "通过" if feo1 > feo0 + 2 else "失败!"))

# 3) 标定链路 (3 炉合成: 输入微扰, 实测取论文值; 仅证明链路可跑)
h1 = base
h2 = copy.deepcopy(base); h2['铁水温度']=1290.0; h2['w([C])']=4.30; h2['终点T_actual']=1660.0; h2['终点C_actual']=0.042
h3 = copy.deepcopy(base); h3['铁水温度']=1260.0; h3['w([Si])']=0.45; h3['终点T_actual']=1648.0; h3['终点C_actual']=0.035
from calibration.optimizer import Calibrator
emit("[4] 标定链路冒烟 (3炉合成数据, maxiter=8, 非真实标定)...")
rep = Calibrator().calibrate([h1,h2,h3], maxiter=8, verbose=False)
emit("    loss=%.3f, 训练集: %s" % (rep['loss'], rep['训练集']))
if '验证集' in rep:
    emit("    验证集: %s" % rep['验证集'])
emit("    标定后参数(节选): eta_post=%.3f, Q_rad_W=%.2e, A_em=%.0f" % (
    rep['参数']['eta_post'], rep['参数']['Q_rad_W'], rep['参数']['A_emulsion']))
emit()
emit("管线冒烟测试完成。真实标定命令: python main.py calibrate data/raw/heats.csv")
emit("(heats.csv 按 data/raw/heats_template.csv 格式, 需 终点C_actual/终点T_actual, 建议含 终点FeO_actual)")

out_path = os.path.join(ROOT, 'data', 'results', 'calib_smoke_output.txt')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(lines) + "\n")
print("\n[结果已保存到 %s]" % out_path)
