import pytest
from model.engine import Engine


@pytest.fixture
def engine():
    return Engine()


@pytest.fixture
def sample_input():
    return {
        '铁水质量': 170.0, '铁水温度': 1276.0,
        'w([C])': 4.48, 'w([Si])': 0.34, 'w([Mn])': 0.19,
        'w([P])': 0.117, 'w([S])': 0.01, 'w([V])': 0.0, 'w([Ti])': 0.0,
        '废钢加入量': 17.0, '供氧流量': 850.0, '底吹氩气流量': 11.5,
        '炉渣碱度': 3.0, 'MgO目标含量': 8.0, '冷却剂类型': 'sinter',
        '终点碳目标': 0.04, '终点温度目标': 1650.0,
        '枪位高度': 1.5, '喷嘴直径': 0.05, '喷嘴夹角': 12.0, '熔池深度': 1.8,
    }


def test_engine_runs(engine, sample_input):
    result = engine.run(sample_input)
    assert '终点C' in result
    assert '终点T' in result
    assert result['终点C'] > 0
    assert result['终点T'] > 1400


def test_engine_produces_time_series(engine, sample_input):
    result = engine.run(sample_input)
    assert '时序数据' in result
    ts = result['时序数据']
    assert len(ts['time']) > 0
    assert len(ts['T']) > 0


def test_endpoint_C_approaches_target(engine, sample_input):
    result = engine.run(sample_input)
    assert result['终点C'] <= sample_input['终点碳目标'] + 0.5
