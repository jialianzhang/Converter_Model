"""输出结果"""
import pandas as pd

class OutputWriter:
    @staticmethod
    def write_results(filepath, results):
        df = pd.DataFrame(results)
        df.to_csv(filepath, index=False, encoding='utf-8-sig')

    @staticmethod
    def write_time_series(filepath, ts_data):
        df = pd.DataFrame(ts_data)
        df.to_csv(filepath, index=False, encoding='utf-8-sig')
