"""临床—正式调度接口验证；所有记录合成，绝不加入真实发现/验证队列。"""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from functions import clinical, pipeline
from functions.configuration import load_config, validate_contract
from functions.storage import digest, fingerprint, read_json, save_json


@pytest.fixture
def config():
    return load_config(Path(__file__).resolve().parents[1] / "config.yaml")


def minimal_config(value):
    result = deepcopy(value)
    result["patterns"]["enabled"] = False
    for family in result["statistics"]["families"].values():
        family["enabled"] = False
    result["statistics"]["families"]["nihss"].update(enabled=True, outcomes=["nihss_1a"])
    for cohort in ["discovery", "validation"]:
        result["clinical"]["cohorts"][cohort].update(confirmed=True,
          rules_confirmed={"id_mapping": True, "mrs_death": True}, id={"field":"id","kind":"string"})
    return result


def test_family_targets_and_covariates(config):
    outcomes = pipeline.clinical_outcomes(config)
    assert len(outcomes) == 29
    families = pd.Series([value["family"] for value in outcomes]).value_counts().to_dict()
    assert families == {"nihss":19,"mrs_single":7,"moca":3}
    for value in outcomes:
        if value["name"].startswith("mrs_"):
            assert value["covariates"] == ["age","sex","log_lesion_volume","prestroke_mrs","acute_nihss"]
        elif value["name"].startswith("nihss_"):
            assert value["covariates"] == ["age","sex","log_lesion_volume"]
        else:
            assert value["covariates"] == ["age","sex","log_lesion_volume","education"]
    assert set(config["patterns"]["targets"]) == {"mrs_3m"} | {value["name"] for value in outcomes if value["family"] == "nihss"}
    config["statistics"]["families"]["mrs_longitudinal"]["enabled"] = True
    long = next(value for value in pipeline.clinical_outcomes(config) if value["family"] == "mrs_longitudinal")
    assert long["visit_times"] == [3,6,12,24,36,48,60]
    assert long["type"] == "longitudinal"


def test_fixed_method_flags_are_not_silently_ignored(config):
    config["patterns"]["imputation_uses_outcome"] = True
    with pytest.raises(ValueError, match="固定方法属性"):
        validate_contract(config)
    config["patterns"]["imputation_uses_outcome"] = False
    config["imaging"]["healthy_contract"]["nodes"] = 134
    with pytest.raises(ValueError, match="固定方法属性"):
        validate_contract(config)


def test_stage_parameters_follow_clinical_requirements(config):
    baseline_image = fingerprint(pipeline.stage_parameters("lesion", config))
    baseline_clinical = fingerprint(pipeline.stage_parameters("clinical", config))
    config["statistics"]["mi"]["m"] = 50
    assert fingerprint(pipeline.stage_parameters("lesion", config)) == baseline_image
    assert fingerprint(pipeline.stage_parameters("clinical", config)) == baseline_clinical
    config["statistics"]["families"]["moca"]["enabled"] = False
    assert fingerprint(pipeline.stage_parameters("clinical", config)) != baseline_clinical


def test_prepared_table_and_statistics_job(tmp_path, config, monkeypatch):
    config = minimal_config(config)
    (tmp_path / "clinical/discovery").mkdir(parents=True)
    (tmp_path / "clinical/validation").mkdir(parents=True)
    (tmp_path / "reference").mkdir()
    save_json(tmp_path / "reference/manifest.json", {"reference_id":"synthetic_reference"})
    sources = {"discovery":pd.DataFrame({"id":["0001","shared"],"AGE":[50,65],"GENDER":[1,2],"A_NIH1A":[0,1]}),
               "validation":pd.DataFrame({"id":["shared","v003"],"AGE":[65,70],"GENDER":[2,1],"A_NIH1A":[1,2]})}
    for name in sources:
        (tmp_path / f"clinical/{name}/{name}.sas7bdat").write_bytes(name.encode())
    monkeypatch.setattr(clinical, "read_sas", lambda path:(sources[path.parent.name].copy(), {}))
    feature_dir = tmp_path / "lesion_stage"
    feature_dir.mkdir()
    images = pd.DataFrame({"subject_id":["0001","v003"], **{f"{channel}_{mode:03d}":[.02,.03] for channel in ["g","w"] for mode in range(1,130)}})
    images.to_csv(feature_dir / "mode_features.tsv.gz",sep="\t",index=False)
    pd.DataFrame({"subject_id":images.subject_id,"lesion_volume_mm3":[100,200]}).to_csv(feature_dir / "summary.tsv",sep="\t",index=False)
    output = tmp_path / "clinical_stage"
    discovery = clinical.prepare_clinical(tmp_path, config, feature_dir, output / "discovery", "discovery")
    validation = clinical.prepare_clinical(tmp_path, config, feature_dir, output / "validation", "validation")
    assert discovery["subject_ids"] == ["0001"]
    assert validation["subject_ids"] == ["v003"]
    # shared 无有效影像，仍不能绕过原始队列成员的外部验证独立性检查。
    with pytest.raises(ValueError, match="重叠"):
        clinical.assert_disjoint(set(discovery["cohort_subject_ids"]), set(validation["cohort_subject_ids"]))
    job = pipeline.statistics_job(tmp_path,"imputation",config,{"clinical":str(output)},tmp_path/"mi_stage")
    assert Path(job["input_csv"]).is_file()
    assert job["reference_id"] == "synthetic_reference"
    assert job["feature_columns"] == discovery["feature_columns"]
    assert job["outcomes"][0]["covariates"] == ["age","sex","log_lesion_volume"]
    imported = pd.read_csv(job["input_csv"], sep="\t",dtype={"subject_id":str})
    assert imported.subject_id.iloc[0] == "0001"
    assert imported.log_lesion_volume.iloc[0] == pytest.approx(np.log1p(100))


def test_id_and_member_files_are_fingerprinted(tmp_path, config):
    config = minimal_config(config)
    clinical_dir = tmp_path / "clinical"
    (clinical_dir / "discovery").mkdir(parents=True)
    (clinical_dir / "validation").mkdir(parents=True)
    for name in ["discovery.sas7bdat","mapping.tsv","members.tsv"]:
        (clinical_dir / "discovery" / name).write_text("before",encoding="utf-8")
    config["clinical"]["cohorts"]["discovery"]["member_file"] = "members.tsv"
    config["clinical"]["cohorts"]["discovery"]["mapping_file"] = "mapping.tsv"
    before = pipeline.stage_sources(tmp_path,"clinical",config,SimpleNamespace(subjects=None))
    assert len(before) == 3
    (clinical_dir / "discovery" / "mapping.tsv").write_text("after",encoding="utf-8")
    after = pipeline.stage_sources(tmp_path,"clinical",config,SimpleNamespace(subjects=None))
    assert fingerprint(before) != fingerprint(after)


def test_python_adapter_changes_invalidate_models_not_lesions(tmp_path, monkeypatch):
    """Python 定义的结局/协变量变化不能复用旧 R 模型；影像算法不受影响。"""
    functions = tmp_path / 'functions'
    functions.mkdir()
    (tmp_path / 'statistics').mkdir()
    for name in ['cohort','space','patient','index','spectral','reference',
                 'clinical','configuration','annotation','spatial_export','plotting','pipeline','results_report']:
        (functions / f'{name}.py').write_text('before',encoding='utf-8')
    (tmp_path / 'stage.py').write_text('stage',encoding='utf-8')
    (tmp_path / 'statistics/run.R').write_text('r code',encoding='utf-8')
    monkeypatch.setattr(pipeline,'r_versions',lambda root:'fixed R packages')
    stages = ['lesion','clinical','imputation','association','patterns','validation','annotation']
    before = {name:pipeline.code_signature(tmp_path,name) for name in stages}
    (functions / 'pipeline.py').write_text('changed outcome and covariate adapter',encoding='utf-8')
    after = {name:pipeline.code_signature(tmp_path,name) for name in stages}
    assert before['lesion'] == after['lesion']
    assert all(before[name] != after[name] for name in stages if name != 'lesion')


def seeded_pipeline(tmp_path, config, monkeypatch):
    """构造带真实校验值的微型归档；不执行影像计算或 R 模型。"""
    monkeypatch.setattr(pipeline,'stage_sources',lambda root,stage,config,args:{'source':stage})
    monkeypatch.setattr(pipeline,'code_signature',lambda root,stage:{'code':stage})
    monkeypatch.setattr(pipeline.signal,'signal',lambda *values:None)
    state = {'status':'completed','stages':{}}
    directory = tmp_path / config['paths']['results'] / 'regression'
    for stage in ['lesion','clinical','imputation','association']:
        upstream = {name:state['stages'][name]['signature'] for name in pipeline.dependencies[stage]}
        if stage == 'association': upstream['clinical_path'] = state['stages']['clinical']['signature']
        signature = fingerprint({'parameters':pipeline.stage_parameters(stage,config),
            'sources':{'source':stage},'code':{'code':stage},'upstream':upstream})
        target = directory / 'stages' / stage / signature[:16]
        target.mkdir(parents=True)
        output = target / ('associations.csv' if stage == 'association' else 'data.txt')
        output.write_text('outcome,mode,q,family,status\nold_unique_result,1,0.02,nihss,estimated\n'
                          if stage == 'association' else stage,encoding='utf-8')
        save_json(target / 'complete.json',{'signature':signature,'outputs':{output.name:digest(output)}})
        state['stages'][stage] = {'status':'completed','signature':signature,'path':str(target)}
    save_json(directory / 'stages.json',state)
    args = SimpleNamespace(run='regression',stages='report',resume=True,subjects=None,token='test-token')
    return directory,state,args


def test_report_only_rechecks_changed_statistics_and_changes_signature(tmp_path, config, monkeypatch):
    directory,state,args = seeded_pipeline(tmp_path,config,monkeypatch)
    assert pipeline.run_pipeline(tmp_path,args,config) == 0
    before = read_json(directory / 'stages.json')['stages']['report']['signature']
    assert 'old_unique_result' in (directory / 'report.html').read_text(encoding='utf-8')
    config['statistics']['mi']['m'] += 1
    assert pipeline.run_pipeline(tmp_path,args,config) == 0
    after = read_json(directory / 'stages.json')
    assert after['stages']['lesion']['status'] == 'completed'
    assert after['stages']['clinical']['status'] == 'completed'
    assert after['stages']['imputation']['status'] == 'stale'
    assert after['stages']['association']['status'] == 'stale'
    assert after['stages']['report']['signature'] != before
    assert 'old_unique_result' not in (directory / 'report.html').read_text(encoding='utf-8')


def test_modified_upstream_product_invalidates_downstream_and_report(tmp_path, config, monkeypatch):
    directory,state,args = seeded_pipeline(tmp_path,config,monkeypatch)
    target = Path(state['stages']['clinical']['path']) / 'data.txt'
    target.write_text('modified clinical table',encoding='utf-8')
    with pytest.raises(ValueError,match='校验改变'):
        pipeline.ensure_stage_fresh(tmp_path,'association',config,args,state)
    assert pipeline.run_pipeline(tmp_path,args,config) == 0
    after = read_json(directory / 'stages.json')
    for stage in ['clinical','imputation','association']:
        assert after['stages'][stage]['status'] == 'stale'
    assert 'old_unique_result' not in (directory / 'report.html').read_text(encoding='utf-8')


def test_resume_publishes_current_progress_before_expensive_source_checks(tmp_path, config, monkeypatch):
    """后台握手后即应读到本次状态，输入哈希检查不能暴露旧失败记录。"""
    directory, state, args = seeded_pipeline(tmp_path, config, monkeypatch)
    args.stages = "lesion"
    state["stages"]["lesion"].update(status="failed", reason="OLD_SCOPE_ERROR")
    save_json(directory / "stages.json", state)
    save_json(directory / "progress.json", {"status": "incomplete", "heartbeat": 1, "token": "old"})
    original_sources = pipeline.stage_sources
    checks = []

    def sources(root, stage, config, arguments):
        progress = read_json(directory / "progress.json")
        control = read_json(directory / "control.json")
        assert progress["status"] == "running"
        assert progress["token"] == control["token"] == args.token
        assert progress["heartbeat"] >= control["started"]
        if not checks:
            assert progress["stages"]["lesion"]["status"] == "pending"
            assert "reason" not in progress["stages"]["lesion"]
        checks.append(stage)
        return original_sources(root, stage, config, arguments)

    monkeypatch.setattr(pipeline, "stage_sources", sources)
    assert pipeline.run_pipeline(tmp_path, args, config) == 0
    assert checks


def test_immediate_stop_tolerates_process_exit_race(tmp_path, monkeypatch):
    """poll 后、killpg 前退出的自有子进程仍按停止处理。"""
    class process:
        pid = 99999
        returncode = 0
        def poll(self): return None
        def wait(self,timeout=None): return 0
    monkeypatch.setattr(pipeline.subprocess,'Popen',lambda *args,**kwargs:process())
    def already_exited(*args): raise ProcessLookupError('already exited')
    monkeypatch.setattr(pipeline.os,'killpg',already_exited)
    result = pipeline.run_child(['fake'],tmp_path,{'graceful':True,'immediate':True},
                                tmp_path/'missing_stop.json','token',lambda:None)
    assert result == 130


def test_preflight_added_inputs_change_signature_without_config_change(tmp_path, config, monkeypatch):
    """预检缺失输入可正常记账；补齐 SAS/ID 映射/病灶后不能复用旧缺失报告。"""
    from functions import cohort
    config = minimal_config(config)
    clinical_dir = tmp_path / 'clinical'
    (clinical_dir / 'discovery').mkdir(parents=True)
    (clinical_dir / 'validation').mkdir(parents=True)
    (tmp_path / 'reference').mkdir()
    config['clinical']['cohorts']['discovery']['member_file'] = 'members.tsv'
    config['clinical']['cohorts']['discovery']['mapping_file'] = 'mapping.tsv'
    def discover(root,selected,options,cohort='discovery'):
        path = root / 'lesion.nii.gz'
        if not path.is_file(): raise ValueError('没有匹配的病灶')
        return [{'subject_id':('0001' if cohort=='discovery' else '0002'),'path':str(path),'sha256':digest(path)}]
    monkeypatch.setattr(cohort,'discover',discover)
    args = SimpleNamespace(subjects=None)
    current = pipeline.stage_sources(tmp_path,'preflight',config,args)
    assert current['lesion']['status'] == 'unavailable'
    assert all(not value['exists'] for value in current['clinical'].values())
    for path in [clinical_dir/'discovery'/'discovery.sas7bdat',clinical_dir/'validation'/'validation.sas7bdat',
                 clinical_dir/'discovery'/'mapping.tsv',clinical_dir/'discovery'/'members.tsv',
                 tmp_path/'lesion.nii.gz',tmp_path/'reference/manifest.json']:
        path.write_bytes(b'engineering test input')
        added = pipeline.stage_sources(tmp_path,'preflight',config,args)
        assert fingerprint(added) != fingerprint(current)
        current = added
    assert current['lesion']['status'] == 'discovered'
    (clinical_dir/'discovery'/'discovery.sas7bdat').write_bytes(b'changed input content')
    assert fingerprint(pipeline.stage_sources(tmp_path,'preflight',config,args)) != fingerprint(current)


@pytest.mark.parametrize("subjects", [None, "all"])
@pytest.mark.parametrize("stages,cached_preflight", [
    ("preflight,lesion", False),
    ("preflight,lesion", True),
    ("lesion", False),
])
def test_lesion_dispatch_with_cached_or_omitted_preflight(
    tmp_path, config, monkeypatch, subjects, stages, cached_preflight
):
    """复用预检或直接运行病灶阶段时，仍须派发两个队列。"""
    from functions import reference

    monkeypatch.setattr(pipeline, "stage_sources", lambda root, stage, config, args: {"source": stage})
    monkeypatch.setattr(pipeline, "code_signature", lambda root, stage: {"code": stage})
    monkeypatch.setattr(pipeline.signal, "signal", lambda *values: None)
    monkeypatch.setattr(reference, "verify_reference", lambda root: {})
    monkeypatch.setattr(reference, "load_reference", lambda root, with_index: {})
    monkeypatch.setattr(reference, "numerical_checks", lambda value: {"synthetic": True})
    checked_cohorts = []

    def preflight_clinical(root, config, cohort):
        checked_cohorts.append(cohort)
        return []

    monkeypatch.setattr(clinical, "preflight_clinical", preflight_clinical)
    dispatched = []

    def run_child(command, directory, *args):
        job = read_json(directory / "job_config.json")
        assert job["active_cohort"] == directory.name
        dispatched.append(directory.name)
        (directory / "synthetic_result.txt").write_text(directory.name, encoding="utf-8")
        save_json(directory / "progress.json", {"status": "completed"})
        return 0

    monkeypatch.setattr(pipeline, "run_child", run_child)
    args = SimpleNamespace(run="dispatch_regression", stages="preflight", resume=False,
                           subjects=subjects, token="test-token")
    directory = tmp_path / config["paths"]["results"] / args.run
    if cached_preflight:
        assert pipeline.run_pipeline(tmp_path, args, config) == 0
        before = read_json(directory / "stages.json")["stages"]["preflight"]
        completion = Path(before["path"]) / "complete.json"
        cached_output = (completion.read_bytes(), completion.stat().st_mtime_ns)

    args.stages = stages
    args.resume = cached_preflight
    assert pipeline.run_pipeline(tmp_path, args, config) == 0
    state = read_json(directory / "stages.json")
    assert state["status"] == "completed"
    assert state["stages"]["lesion"]["status"] == "completed"
    assert dispatched == ["discovery", "validation"]
    assert checked_cohorts == (["discovery", "validation"] if "preflight" in stages else [])
    if cached_preflight:
        assert cached_output == (completion.read_bytes(), completion.stat().st_mtime_ns)

    # 已完成病灶再次续跑须复用产物，不重新派发子进程。
    args.resume = True
    assert pipeline.run_pipeline(tmp_path, args, config) == 0
    assert dispatched == ["discovery", "validation"]
