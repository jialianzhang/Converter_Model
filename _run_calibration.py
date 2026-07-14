# -*- coding: utf-8 -*-
"""100炉数据检查 + 参数标定 (V5 第五阶段)

用法:
    python _run_calibration.py check            # 仅数据检查, 报告存 data/results/data_check_report.txt
    python _run_calibration.py                  # 检查 + 标定 (maxiter=250, 并行=CPU-2)
    python _run_calibration.py 150 4            # 自定义: maxiter=150, 并行进程=4

标定结果存 data/calibrated/calibration_result.txt。预计用时: 并行6进程约20~35分钟。
终止模式: 氧锚定 (吹入每炉实际耗氧量后读终点, 避免"碳目标=实测碳"的循环验证)。
"""
import sys, io, os, time
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from iomodules.input_reader import InputReader

CSV = os.path.join(ROOT, 'data', 'raw', 'heats_100.csv')

# (键, 下限, 上限) —— 超范围仅警告并列出炉号, 不自动剔除
RANGE_CHECKS = [
    ('铁水质量', 100, 200), ('铁水温度', 1150, 1450),
    ('w([C])', 3.0, 5.5), ('w([Si])', 0.10, 1.0), ('w([Mn])', 0.05, 0.5),
    ('w([P])', 0.05, 0.25), ('w([S])', 0.005, 0.08), ('废钢加入量', 10, 60),
    ('供氧流量', 400, 900), ('终点C_actual', 0.02, 0.12), ('终点T_actual', 1550, 1720),
    ('终点O_actual_ppm', 200, 1500), ('出钢量_actual', 130, 190),
    ('实际石灰kg', 2000, 9000), ('实际白云石kg', 500, 5000), ('实际耗氧量Nm3', 6000, 9500),
]


def load_and_check():
    lines = []
    def emit(s=""):
        print(s); lines.append(s)

    raw = InputReader.read_csv(CSV)
    heats = [InputReader.parse_single_heat(r) for r in raw]
    ids = [str(r.get('炉号', i + 1)) for i, r in enumerate(raw)]
    emit(f"数据文件: {CSV}")
    emit(f"炉次数: {len(heats)}")
    ok = True

    # 1) 标定必需字段
    miss = [ids[i] for i, h in enumerate(heats)
            if h['终点C_actual'] is None or h['终点T_actual'] is None or h['实际耗氧量Nm3'] is None]
    if miss:
        emit(f"[不通过] 缺终点C/T实测或实际耗氧量的炉次: {miss}")
        ok = False
    else:
        emit("[通过] 终点C_actual/终点T_actual/实际耗氧量Nm3 全部齐备 (氧锚定终止可用)")

    # 2) 范围检查
    n_warn = 0
    for key, lo, hi in RANGE_CHECKS:
        vals = [(ids[i], h.get(key)) for i, h in enumerate(heats) if h.get(key) is not None]
        if not vals:
            emit(f"[提示] 列缺失: {key} (将用默认值/替代方案)")
            continue
        bad = [(n, v) for n, v in vals if not (lo <= v <= hi)]
        arr = [v for _, v in vals]
        stat = f"min={min(arr):.4g} max={max(arr):.4g} mean={sum(arr)/len(arr):.4g}"
        if bad:
            n_warn += len(bad)
            emit(f"[警告] {key} 超范围[{lo},{hi}] {len(bad)}炉: {bad[:5]}{'...' if len(bad)>5 else ''} | {stat}")
        else:
            emit(f"[通过] {key}: {stat}")

    # 3) 交叉一致性
    bad_t, bad_w, bad_co = [], [], []
    for i, h in enumerate(heats):
        if h.get('吹炼时长s') is not None:
            t_calc = h['实际耗氧量Nm3'] / h['供氧流量'] * 60.0
            if abs(t_calc - h['吹炼时长s']) > 90:
                bad_t.append(ids[i])
        if h.get('出钢量_actual') is not None:
            ratio = h['出钢量_actual'] / (h['铁水质量'] + h['废钢加入量'])
            if not (0.84 <= ratio <= 0.98):
                bad_w.append((ids[i], round(ratio, 3)))
        if h.get('终点O_actual_ppm') is not None:
            prod = h['终点C_actual'] * h['终点O_actual_ppm'] / 1e4
            if not (0.0012 <= prod <= 0.0070):
                bad_co.append((ids[i], round(prod, 5)))
    emit(f"[{'通过' if not bad_t else '警告'}] 耗氧量/流量 vs 吹炼时长 一致性: 异常 {bad_t or '无'}")
    emit(f"[{'通过' if not bad_w else '警告'}] 出钢量/(铁水+废钢) ∈[0.84,0.98]: 异常 {bad_w or '无'}")
    emit(f"[{'通过' if not bad_co else '警告'}] 终点C×O 碳氧积 ∈[0.0012,0.0070]: 异常 {bad_co or '无'}")

    # 4) 数据集级假设与默认值声明
    emit()
    emit("数据集级处理声明:")
    emit("  - 终点碳目标列 = 终点C_actual (循环) → 已改用氧锚定终止, 该列不参与标定")
    emit("  - 实际烧结矿kg 缺失 → 按 0 处理 (大废钢比操作, 假设无矿石冷却剂; 若实际有加请修正)")
    emit("  - 底吹氩缺失 → 厂定额默认 11.5 m³/min; 枪位/喷嘴/熔池深度 → 默认常数 (ε比值中抵消)")
    emit("  - w_V/w_Ti 缺失 → 0; 碱度/MgO目标 → 不被消费 (实际辅料量已提供)")
    for h in heats:
        if h.get('实际烧结矿kg') is None:
            h['实际烧结矿kg'] = 0.0

    emit()
    emit(f"结论: {'数据满足标定要求' if ok else '存在阻断问题, 修正后重跑'} (警告 {n_warn} 项为提示性, 不阻断)")

    out = os.path.join(ROOT, 'data', 'results', 'data_check_report.txt')
    with open(out, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n[检查报告已保存到 {out}]")
    return heats, ok


if __name__ == '__main__':
    heats, ok = load_and_check()
    if len(sys.argv) > 1 and sys.argv[1] == 'check':
        sys.exit(0 if ok else 1)
    if not ok:
        print("数据检查未通过, 终止。")
        sys.exit(1)

    maxiter = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    try:
        import multiprocessing
        default_jobs = max(1, multiprocessing.cpu_count() - 2)
    except Exception:
        default_jobs = 1
    n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else default_jobs

    print(f"\n=== 参数标定: {len(heats)} 炉, maxiter={maxiter}, 并行={n_jobs} ===")
    print("(氧锚定终止; 指标: C/T 必选 + 终点氧a_O(w=0.5) + 出钢量(w=0.5))")
    t0 = time.time()
    from calibration.optimizer import Calibrator
    report = Calibrator().calibrate(heats, maxiter=maxiter, n_jobs=n_jobs)
    elapsed = time.time() - t0

    lines = [f"=== 标定结果 (用时 {elapsed/60:.1f} min) ==="]
    for k, v in report.items():
        lines.append(f"{k}: {v}")
    text = "\n".join(lines)
    print(text)
    out = os.path.join(ROOT, 'data', 'calibrated', 'calibration_result.txt')
    with open(out, 'w', encoding='utf-8') as f:
        f.write(text + "\n")
    print(f"\n[标定结果已保存到 {out}]")
