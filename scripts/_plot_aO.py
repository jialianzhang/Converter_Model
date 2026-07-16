# -*- coding: utf-8 -*-
import sys,io,os; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,ROOT)   # scripts/ 上一级 = 项目根
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')
import numpy as np
from model.engine import Engine
from thermo.wagner import calc_lgf_i
from iomodules.input_reader import InputReader
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def curve(h,eta=0.45):
    h2=dict(h); h2['终止模式']='目标碳'; h2['终点碳目标']=h.get('终点C_actual',0.04)
    if h2.get('实际烧结矿kg') is None: h2['实际烧结矿kg']=0.0
    ts=Engine(params={'eta_O2_late':eta}).run(h2)['时序数据']
    t=np.array(ts['time_s']);aO=np.array(ts['a_O_star']);TK=np.array(ts['T_K']);wc=np.array(ts['w_C'])
    ppm=np.array([aO[i]/max(10.0**calc_lgf_i('O',{'C':ts['w_C'][i],'Si':ts['w_Si'][i],'Mn':ts['w_Mn'][i],'P':ts['w_P'][i],'S':ts['w_S'][i]},TK[i]),1e-6)*1e4 for i in range(len(t))])
    coeq=np.array([0.0025/max(wc[i],1e-4)*1e4 for i in range(len(t))])  # C-O平衡本体氧ppm
    return t,aO,ppm,wc,coeq
z={'铁水质量':170.0,'铁水温度':1276.0,'w([C])':4.48,'w([Si])':0.34,'w([Mn])':0.19,'w([P])':0.117,'w([S])':0.01,'w([V])':0.0,'w([Ti])':0.0,'废钢加入量':17.0,'供氧流量':850.0,'底吹氩气流量':11.5,'炉渣碱度':3.0,'MgO目标含量':8.0,'终点碳目标':0.04,'终点温度目标':1650.0,'枪位高度':1.5,'喷嘴直径':0.05,'喷嘴夹角':12.0,'熔池深度':1.8,'终点C_actual':0.04}
t,aO,ppm,wc,coeq=curve(z)
heats=[InputReader.parse_single_heat(r) for r in InputReader.read_csv(os.path.join(ROOT,'data','raw','heats_100.csv'))]
mp=[];ap=[]
for h in heats:
    if h.get('终点O_actual_ppm') is None: continue
    _,_,pp,_,_=curve(h); mp.append(pp[-1]); ap.append(h['终点O_actual_ppm'])
mp=np.array(mp);ap=np.array(ap)
fig,ax=plt.subplots(1,3,figsize=(16,4.6))
a0=ax[0];a0b=a0.twinx()
a0.plot(t/60,aO,'b-',lw=2,label='a_O* (activity)')
a0.set_yscale('log'); a0.set_xlabel('time (min)'); a0.set_ylabel('a_O* activity',color='b'); a0.grid(alpha=0.3)
a0b.plot(t/60,wc,'g--',lw=1.3,label='w_C %'); a0b.set_ylabel('w_C (%)',color='g')
a0.set_title('Zhang heat: interface a_O* vs time (endgame explosion)')
ax[1].plot(t/60,ppm,'r-',lw=2,label='model interface [O]')
ax[1].plot(t/60,coeq,'k--',lw=1.5,label='C-O equilib bulk [O] ([C][O]=0.0025)')
ax[1].set_xlabel('time (min)'); ax[1].set_ylabel('[O] ppm'); ax[1].set_ylim(0,3000); ax[1].grid(alpha=0.3)
ax[1].legend(); ax[1].set_title('interface [O] ppm vs C-O bulk equilibrium')
ax[2].scatter(ap,mp,s=20,alpha=0.6); lim=max(mp.max(),ap.max())*1.05
ax[2].plot([0,lim],[0,lim],'k--',alpha=0.6,label='y=x (model=actual)')
ax[2].set_xlabel('actual endpoint [O] TSO bulk (ppm)'); ax[2].set_ylabel('model endpoint interface [O] (ppm)')
ax[2].legend(); ax[2].grid(alpha=0.3); ax[2].set_title('100 heats: model interface vs actual bulk O (1.5x)')
plt.tight_layout(); out=os.path.join(ROOT,'data','results','a_O_curve.png'); plt.savefig(out,dpi=110)
print("[图已重画(英文标签)]",out)
