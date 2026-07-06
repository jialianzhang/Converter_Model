# 转炉冶炼终点预测模型 V2 — 计算流程与核心公式

> 面向PPT展示，按计算流程组织，每步仅列出关键公式。

---

## 整体计算流程

```
┌──────────────┐    ┌─────────────────────────────────────────────────────────────┐
│  ① 静态模型   │    │  ② 时间步循环 (13-step Euler-forward, Δt自适应1~5s)          │
│  辅料·氧耗·   │    │                                                             │
│  出钢量计算    │    │  Step 1-2    Step 3       Step 4-6      Step 7-8           │
│              │    │  混合能+传质 → 平衡常数+活度 → 耦合反应求解 → 物质平衡更新    │
│              │    │                                                             │
│              │    │  Step 9-10    Step 11-12    终止判断                         │
│              │    │  辅料+废钢熔化 → 热量平衡+温度更新 → C≤目标? → 输出           │
└──────────────┘    └─────────────────────────────────────────────────────────────┘
```

---

## 一、静态模型 — 初始条件计算

| 步骤 | 计算内容 | 核心公式 |
|------|---------|---------|
| 物料平衡 | 元素氧化量→渣量→辅料需求 | $\sum m_{in} = \sum m_{out}$ |
| 石灰加入量 | 碱度约束 CaO/SiO₂ | $m_{lime} = \dfrac{R \cdot m_{SiO_2}}{w_{CaO}^{lime}}$ |
| 白云石加入量 | MgO目标约束 | $m_{dolo} = f(MgO_{target})$ |
| 耗氧量 | 元素氧化耗O₂ | $V_{O_2} = \sum \dfrac{\Delta m_i}{M_i} \cdot n_i \cdot \dfrac{0.016}{1.429}$ |
| 出钢量+铁耗 | 铁平衡 | $\eta_{Fe} = \dfrac{m_{铁水}+m_{废钢}}{m_{出钢}} \times 1000$ |

---

## 二、Step 1 — 混合能计算 (Nakanishi)

| 计算内容 | 核心公式 |
|---------|---------|
| 顶吹混合能 | $\varepsilon_T = \dfrac{0.0403 \cdot Q_T \cdot D \cdot U^2 \cdot \cos^2\theta}{W \cdot X}$ |
| 底吹混合能 | $\varepsilon_B = \dfrac{28.5 \cdot Q_B \cdot T}{W} \cdot \lg\!\left(1 + \dfrac{H}{148}\right)$ |
| 总混合能 | $\varepsilon_{total} = \varepsilon_T + \varepsilon_B$ |

---

## 三、Step 2 — 动态修正传质系数

| 计算内容 | 核心公式 |
|---------|---------|
| 钢液传质系数 | $k_m = k_m^0 \left( \dfrac{\varepsilon_{total}}{\varepsilon_{total}^0} \right)^{\gamma}$ |
| 炉渣传质系数 | $k_s' = k_s^0 \cdot \min\!\left(1 + \dfrac{3q_{CO}}{q_{Ar}}, 10\right) \cdot \dfrac{40}{m_{slag}}$ |
| 固相率修正 | $k_s = k_s' \times 10^{-n \cdot f_{solid}}$ |
| 钢液侧F系数 | $F_M = \dfrac{k_m \cdot \rho_m}{100 \cdot M_M}$ |
| 炉渣侧F系数 | $F_{MO_n} = \dfrac{k_s \cdot \rho_s}{100 \cdot M_{MO_n}}$ |
| 反应面积 | $A = A_{eff} \times (1 - f_{solid})^{2/3}$ |

---

## 四、Step 3 — 平衡常数与活度

| 计算内容 | 核心公式 |
|---------|---------|
| 反应平衡常数 | $\lg K_M = \dfrac{A}{T} + B$ |
| 表观平衡常数 | $E_M = \dfrac{100 \cdot C \cdot N_{MO_n} \cdot f_M \cdot K_M}{\rho_S \cdot f_{MO_n}}$ |
| 金属活度 (Wagner) | $\lg f_i = \sum_j e_i^j \cdot w([j])$ |
| 炉渣活度 (离子理论) | $f_{iO_n} = \dfrac{1}{\sum_j \frac{w_j}{100} \cdot \exp(-\varepsilon_{ij}/RT)}$ |

---

## 五、Step 4-6 — 耦合反应求解器（核心）

| 计算内容 | 核心公式 |
|---------|---------|
| 传质通量 | $J_M = F_M \{ w([M])^b - w([M])^* \}$ |
| 脱碳通量 | $J_C = F_C \{ w([C])^b - w([C])^* \} = G_{CO}(P_{CO}^* - 1)$ |
| **氧平衡方程** | $2J_{Si} + J_{Mn} + 2.5J_P + 1.5J_V + 2J_{Ti} + J_{Fe} + J_C - J_O = 0$ |
| O₂注入限制 | $J_{O2}^{inj} \propto \max(0.04,\ w_C/0.8)$ |
| 渣侧氧势联动 | $a_O^{cap} = \max(0.3,\ 1.5 \cdot w_C/0.5)$ |
| 界面浓度 | $w([M])^* = \dfrac{F_{MO_n} \cdot w((MO_n))^b + F_M \cdot w([M])^b}{F_M + F_{MO_n} \cdot E_M \cdot (a_O^*)^n}$ |
| 碳专属界面浓度 | $w([C])^* = \dfrac{F_C \cdot w([C])^b + G_{CO}}{F_C + G_{CO} \cdot K_C \cdot a_O^*}$ |
| 反应速率 | $-\dfrac{\Delta m_i}{\Delta t} = A \cdot F_i \cdot M_i \{ w([i])^b - w([i])^* \}$ |

---

## 六、Step 7-8 — 物质平衡更新

| 计算内容 | 核心公式 |
|---------|---------|
| 钢液元素更新 | $m_{elem}^{(l)} = m_{elem}^{(l-1)} - \Delta m_{elem}$ |
| 炉渣氧化物更新 | $\Delta m_{ox} = \Delta m_{elem} \times f_{conv}$ |
| FeO控制 | 源头氧势调控（Step 4），不截断质量流 |

---

## 七、Step 9 — 辅料溶解 (O₂消耗量驱动)

| 计算内容 | 核心公式 |
|---------|---------|
| O₂进度 | $\phi = O_2^{consumed} / O_2^{total}$ |
| 辅料加入比例 | $f_{flux} = \max(0,\ \min(1,\ (\phi - 0.10)/0.80))$ |
| 石灰溶解 (Kitamura) | $\dot{m}_{lime} = n \cdot A \cdot \left\{ \dfrac{k' \cdot \rho_s}{100 \cdot r'} \right\} \cdot \Delta(\%\mathrm{CaO}) \cdot r' \cdot x$ |

---

## 八、Step 10 — 废钢熔化 (抛物线)

| 计算内容 | 核心公式 |
|---------|---------|
| 抛物线熔化速率 | $\dot{m}_{scrap}(t) = m_{scrap}^{total} \cdot \dfrac{6\tau(1-\tau)}{t_{end}},\ \tau = \dfrac{t}{t_{end}}$ |

---

## 九、Step 11-12 — 热量平衡与温度更新

| 计算内容 | 核心公式 |
|---------|---------|
| 热收入 | $Q_{HG} = \sum \dfrac{|\Delta m_i|}{M_i} \cdot \Delta H_i^{oxid} + \eta_{post} \cdot Q_C^{oxid} + Q_{C_2S}$ |
| 废钢吸热 | $Q_{scrap} = \Delta m_{scrap} \times 1.3 \times 10^6$ |
| 辅料分解吸热 | $Q_{flux} = \Delta m_{lime}{\times}1.78{\times}10^6 + \Delta m_{dolo}{\times}3.02{\times}10^6$ |
| CO炉气显热 | $Q_{offgas} = \dfrac{q_{CO}}{0.0224} \cdot 33 \cdot (T-298) \cdot \Delta t$ |
| 烧结冷却(延期) | $Q_{sinter} = \begin{cases} 0, & w_C \geq 2.5\% \\ \dot{m}_{sinter} \times 1.8\times10^6, & w_C < 2.5\% \end{cases}$ |
| **温度更新** | $\Delta T = \dfrac{Q_{HG} - Q_{consumed}}{c_M m_M + c_S m_S}$ |

---

## 十、自适应步长与终止

| 条件 | 步长 | 终止 |
|------|------|------|
| $w([C]) > 0.8\%$ | $\Delta t = 1$s | |
| $0.3\% < w([C]) \leq 0.8\%$ | $\Delta t = 2$s | |
| $w([C]) \leq 0.3\%$ | $\Delta t = 5$s | |
| $O_2 \geq 50\% \land w([C]) \leq target$ | — | **循环结束** |

---

## 公式汇总索引 (按模块)

| 模块 | 关键公式 | 来源 |
|------|---------|------|
| 混合能 | $\varepsilon_T, \varepsilon_B$ (Nakanishi) | 张强(2024) |
| 动态传质 | $k_m, k_s, F_M, F_{MO_n}$ | Kitamura(2009)+V2工业标定 |
| 活度 | Wagner + 离子结构理论 | 张强(2024) |
| 耦合反应 | $J_M$, 氧平衡(式6), 界面浓度(式7) | 张强(2024)/Kitamura(2009) |
| O₂控制 | $J_{O2}^{inj}$限流 + $a_O^{cap}$联动 | **V2核心修复** |
| 物质平衡 | 钢液/炉渣质量更新(式12) | 张强(2024) |
| 辅料熔化 | Kitamura石灰溶解 + O₂进度驱动 | Kitamura(2009)+V2改造 |
| 废钢熔化 | 抛物线分配 | V2经验模型 |
| 热量平衡 | 分源热收支 + 温度更新(式13) | 张强(2024)+V2分源改造 |
| 寻根 | Brent→NR→Bisection | V2升级 |
