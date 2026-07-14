"""读取生产数据"""
import pandas as pd

class InputReader:
    @staticmethod
    def read_csv(filepath):
        df = pd.read_csv(filepath)
        return df.to_dict('records')

    @staticmethod
    def _opt_float(row, key):
        """可选列: 缺失/空/NaN → None"""
        v = row.get(key, None)
        if v is None:
            return None
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return None if f != f else f   # NaN 检查

    @staticmethod
    def parse_single_heat(row):
        out = {
            '铁水质量': float(row.get('铁水质量', 170)),
            '铁水温度': float(row.get('铁水温度', 1276)),
            'w([C])': float(row.get('w_C', 4.48)),
            'w([Si])': float(row.get('w_Si', 0.34)),
            'w([Mn])': float(row.get('w_Mn', 0.19)),
            'w([P])': float(row.get('w_P', 0.117)),
            'w([S])': float(row.get('w_S', 0.01)),
            'w([V])': float(row.get('w_V', 0.0)),
            'w([Ti])': float(row.get('w_Ti', 0.0)),
            '废钢加入量': float(row.get('废钢', 17)),
            '供氧流量': float(row.get('供氧流量', 850)),
            '底吹氩气流量': float(row.get('底吹氩', 11.5)),
            '炉渣碱度': float(row.get('碱度', 3.0)),
            'MgO目标含量': float(row.get('MgO目标', 8.0)),
            '冷却剂类型': 'sinter',
            '终点碳目标': float(row.get('终点碳目标', 0.04)),
            '终点温度目标': float(row.get('终点温度目标', 1650)),
            '枪位高度': float(row.get('枪位高度', 1.5)),
            '喷嘴直径': float(row.get('喷嘴直径', 0.05)),
            '喷嘴夹角': float(row.get('喷嘴夹角', 12.0)),
            '熔池深度': float(row.get('熔池深度', 1.8)),
        }
        # 标定用实测终点列 (C/T 必填才能参与标定, FeO 强烈建议 —— 防FeO与二次燃烧率补偿)
        out['终点C_actual'] = InputReader._opt_float(row, '终点C_actual')
        out['终点T_actual'] = InputReader._opt_float(row, '终点T_actual')
        out['终点FeO_actual'] = InputReader._opt_float(row, '终点FeO_actual')
        # 可选扩展列: 实际辅料量/耗氧量由引擎优先消费 (替代静态估算);
        # 终点O_actual_ppm(副枪TSO氧)/出钢量_actual 为终渣FeO缺失时的替代约束指标
        for key in ['终点P_actual', '终点Mn_actual', '实际石灰kg', '实际白云石kg',
                    '实际烧结矿kg', '实际耗氧量Nm3', '吹炼时长s', '炉龄',
                    '终点O_actual_ppm', '出钢量_actual']:
            out[key] = InputReader._opt_float(row, key)
        return out
