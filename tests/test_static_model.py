import pytest
from model.static_model import StaticModel

@pytest.fixture
def static_model():
    return StaticModel()

@pytest.fixture
def sample_input():
    return {
        '铁水质量': 170.0, '铁水温度': 1276.0,
        'w([C])': 4.48, 'w([Si])': 0.34, 'w([Mn])': 0.19,
        'w([P])': 0.117, 'w([S])': 0.01, 'w([V])': 0.0, 'w([Ti])': 0.0,
        '废钢加入量': 17.0, '供氧流量': 850.0, '底吹氩气流量': 11.5,
        '炉渣碱度': 3.0, 'MgO目标含量': 8.0, '冷却剂类型': 'sinter',
        '终点碳目标': 0.04, '终点温度目标': 1650.0,
    }

def test_static_model_runs(static_model, sample_input):
    result = static_model.calculate(sample_input)
    assert '总耗氧量' in result
    assert '石灰加入量' in result
    assert result['总耗氧量'] > 0

def test_lime_addition_positive(static_model, sample_input):
    result = static_model.calculate(sample_input)
    assert result['石灰加入量'] > 0

def test_iron_loss_calculation(static_model, sample_input):
    result = static_model.calculate(sample_input)
    expected_bl = (sample_input['铁水质量'] + sample_input['废钢加入量']) / result['出钢量'] * 1000
    assert abs(result['铁耗'] - expected_bl) < 1.0

def test_oxygen_consumption_reasonable(static_model, sample_input):
    result = static_model.calculate(sample_input)
    assert 5000 < result['总耗氧量'] < 12000
