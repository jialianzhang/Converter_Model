# -*- coding: utf-8 -*-
"""最终标定 (V1.0 冻结前, 目标碳框架=使用框架)

- framing='target_c': 吹到实测碳停 (张强口径), C 构造性命中 (w_C=0), 温度/Mn/P 为真预测
- eta_O2_late 固定 0.45 (由终渣 FeO 工业范围 15~30% 物理定值, 不参与优化)
- 目标: T(1.0) + Mn(0.5) + P(0.3); 出钢量/aO 仅作报告 (出钢量含未建模的喷溅损失,
  作目标会把铁氧化推向虚高; aO 为界面vs本体偏置)
用法: python _run_final_calibration.py [maxiter=300] [n_jobs=CPU-2]
结果: data/calibrated/final_calibration_v1.txt
"""
import sys, io, os, time
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
from iomodules.input_reader import InputReader

if __name__ == '__main__':
    heats = [InputReader.parse_single_heat(r) for r in InputReader.read_csv(
        os.path.join(ROOT, 'data', 'raw', 'heats_100.csv'))]
    for h in heats:
        if h.get('实际烧结矿kg') is None: h['实际烧结矿kg'] = 0.0

    maxiter = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    try:
        import multiprocessing
        default_jobs = max(1, multiprocessing.cpu_count() - 2)
    except Exception:
        default_jobs = 1
    n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else default_jobs

    print(f"=== 最终标定 (目标碳框架): {len(heats)} 炉, maxiter={maxiter}, 并行={n_jobs} ===")
    print("(eta_O2_late 固定 0.45; 目标 T+Mn+P; C 构造性命中不入目标)")
    t0 = time.time()
    from calibration.optimizer import Calibrator
    report = Calibrator().calibrate(
        heats, maxiter=maxiter, n_jobs=n_jobs,
        framing='target_c', fixed={'eta_O2_late': 0.45},
        w_C=0.0, w_T=1.0, w_FeO=0.0, w_aO=0.0, w_W=0.0, w_Mn=0.5, w_P=0.3)
    elapsed = time.time() - t0

    lines = [f"=== 最终标定结果 V1.0 (用时 {elapsed/60:.1f} min) ==="]
    for k, v in report.items():
        lines.append(f"{k}: {v}")
    text = "\n".join(lines)
    print(text)
    out = os.path.join(ROOT, 'data', 'calibrated', 'final_calibration_v1.txt')
    with open(out, 'w', encoding='utf-8') as f:
        f.write(text + "\n")
    print(f"\n[结果已保存到 {out}]")
