"""读取生产数据"""
import pandas as pd

class InputReader:
    @staticmethod
    def read_csv(filepath):
        df = pd.read_csv(filepath)
        return df.to_dict('records')

    @staticmethod
    def parse_single_heat(row):
        return {
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
