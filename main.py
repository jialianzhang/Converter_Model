"""转炉冶炼终点预测模型 V2 - 主入口"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from model.engine import Engine
from iomodules.input_reader import InputReader
from iomodules.output_writer import OutputWriter

def run_single_heat(input_dict=None):
    """运行单炉次模拟"""
    engine = Engine()
    if input_dict is None:
        input_dict = {
            '铁水质量': 170.0, '铁水温度': 1276.0,
            'w([C])': 4.48, 'w([Si])': 0.34, 'w([Mn])': 0.19,
            'w([P])': 0.117, 'w([S])': 0.01, 'w([V])': 0.0, 'w([Ti])': 0.0,
            '废钢加入量': 17.0, '供氧流量': 850.0, '底吹氩气流量': 11.5,
            '炉渣碱度': 3.0, 'MgO目标含量': 8.0, '冷却剂类型': 'sinter',
            '终点碳目标': 0.04, '终点温度目标': 1650.0,
            '枪位高度': 1.5, '喷嘴直径': 0.05, '喷嘴夹角': 12.0, '熔池深度': 1.8,
        }
    result = engine.run(input_dict)
    print(f"终点C: {result['终点C']:.4f}%")
    print(f"终点T: {result['终点T']:.1f}°C")
    print(f"模拟步数: {result['步数']}")
    print(f"吹炼时间: {result['时序数据']['time_s'][-1]:.0f}s")
    return result

def run_batch(input_csv, output_csv=None):
    """批量运行多炉次"""
    reader = InputReader()
    engine = Engine()
    heats = reader.read_csv(input_csv)
    results = []
    for i, heat in enumerate(heats):
        inputs = reader.parse_single_heat(heat)
        result = engine.run(inputs)
        results.append({
            '炉次': i + 1,
            '终点C_预测': result['终点C'],
            '终点T_预测': result['终点T'],
            '吹炼时间_s': result['时序数据']['time'][-1],
        })
        print(f"炉次 {i+1}/{len(heats)}: C={result['终点C']:.4f}%, T={result['终点T']:.1f}°C")

    if output_csv:
        OutputWriter.write_results(output_csv, results)
        print(f"结果已保存到 {output_csv}")
    return results

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_batch(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        print("=== 转炉冶炼终点预测模型 V2 ===")
        run_single_heat()
