"""报告保留诊断与缺失状态，不读取未完成阶段或测试替代资料。"""
import json
from pathlib import Path
import pandas as pd
from functions.results_report import report_sections, table_section


def test_report_skips_stale_and_preserves_empty(tmp_path):
    empty=tmp_path/'empty.csv'
    empty.touch()
    assert '没有可展示' in table_section(empty,'未估计结果',tmp_path)
    state={'stages':{'association':{'status':'stale','path':str(tmp_path/'absent')}}}
    assert report_sections(tmp_path,state)==''


def test_report_shows_imputation_and_model_diagnostics(tmp_path):
    mi=tmp_path/'mi/mrs_3m'
    mi.mkdir(parents=True)
    (mi/'diagnostics.json').write_text(json.dumps({'status':'mi_warning','m':30,'iterations':100}))
    model=tmp_path/'association'
    model.mkdir()
    pd.DataFrame([{'family':'mrs_single','status':'assumption_limited','p_joint':.2}]).to_csv(model/'associations.csv',index=False)
    pd.DataFrame([{'outcome':'mrs_3m','po_status':'violation','po_min_p':.001}]).to_csv(model/'model_diagnostics.csv',index=False)
    state={'stages':{'imputation':{'status':'completed','path':str(mi.parent)},
                     'association':{'status':'completed','path':str(model)}}}
    result=report_sections(tmp_path,state)
    assert 'mi_warning' in result and 'violation' in result
    assert 'model_diagnostics.csv' in result
