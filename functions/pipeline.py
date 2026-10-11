"""正式阶段调度：按依赖冻结指纹；下游参数改变时保留已完成的病灶计算。"""
from __future__ import annotations
import html
import os
import signal
import subprocess
import sys
import time
import uuid
from functools import lru_cache
import importlib.metadata
from pathlib import Path
from functions.storage import digest, fingerprint, read_json, save_json, verify_files
from functions.control import run_directory, run_lock, identity


stage_order = ['preflight','lesion','clinical','imputation','association','patterns','validation','annotation','report']
dependencies = {'preflight':[], 'lesion':[], 'clinical':['lesion'], 'imputation':['clinical'],
                'association':['imputation'], 'patterns':['clinical'],
                'validation':['patterns','clinical'], 'annotation':['validation','lesion','patterns'], 'report':[]}


def stage_parameters(stage: str, config: dict) -> dict:
    """仅签名本阶段相关设置；临床配置不进入 lesion 的签名。"""
    sections = {'preflight':['paths','subjects','imaging','clinical'],
                'lesion':['imaging','subjects'], 'clinical':['clinical'],
                'imputation':['statistics','clinical'], 'association':['statistics'],
                'patterns':['patterns','statistics','clinical'], 'validation':['validation','statistics'],
                'annotation':['annotation'], 'report':[]}
    parameters = {key:config.get(key) for key in sections[stage]}
    if stage in ['preflight','clinical']:
        parameters['requirements'] = {'families':config['statistics']['families'],
            'patterns':{key:config['patterns'].get(key) for key in ['enabled','targets']}}
    return parameters


def stage_sources(root: Path, stage: str, config: dict, args) -> dict:
    """记录真实输入内容；不把私人数据复制进配置或离线包。"""
    if stage == 'preflight':
        from functions.clinical import clinical_input_files
        from functions.cohort import discover, configured_cohorts
        def file_state(path):
            path = Path(path)
            return {'exists':path.is_file(), 'sha256':digest(path) if path.is_file() else None}
        # 尚未到位的数据也是明确的输入状态，补齐后必须重做预检。
        sources = {'clinical':{str(path):file_state(path) for path in clinical_input_files(root,config)},
                   'healthy_manifest':file_state(root/'reference/manifest.json')}
        selected = args.subjects or config['subjects'].get('default','all')
        if isinstance(selected,list): selected = ','.join(selected)
        try:
            lesion_cohorts = {}
            for cohort in configured_cohorts(config):
                lesion_cohorts[cohort] = discover(root, selected, config, cohort=cohort)
            if "validation" in lesion_cohorts:
                from functions.clinical import assert_disjoint
                assert_disjoint(
                    {item["subject_id"] for item in lesion_cohorts["discovery"]},
                    {item["subject_id"] for item in lesion_cohorts["validation"]},
                )
            sources["lesion"] = {"status": "discovered", "cohorts": lesion_cohorts}
        except (ValueError,FileNotFoundError,OSError) as exc:
            sources["lesion"] = {"status":"unavailable","reason":f"{type(exc).__name__}: {exc}"}
        return sources
    if stage == 'lesion':
        from functions.cohort import discover, configured_cohorts, method_identity
        from functions.reference import verify_reference
        selected = args.subjects or config['subjects'].get('default','all')
        if isinstance(selected,list): selected = ','.join(selected)
        cohorts = {}
        for cohort in configured_cohorts(config):
            cohorts[cohort] = discover(root, selected, config, cohort=cohort)
        if "validation" in cohorts:
            from functions.clinical import assert_disjoint
            assert_disjoint(
                {item["subject_id"] for item in cohorts["discovery"]},
                {item["subject_id"] for item in cohorts["validation"]},
            )
        return {"cohorts": cohorts,
                "method":method_identity(root,config,verify_reference(root))}
    if stage == 'clinical':
        from functions.clinical import clinical_input_files
        return {str(path):digest(path) for path in clinical_input_files(root,config) if Path(path).is_file()}
    if stage == 'annotation':
        path=root/config['annotation']['resources']/'manifest.json'
        return {'annotation_manifest':digest(path)}
    return {}


@lru_cache(maxsize=1)
def r_versions(root: Path) -> str:
    code='for(p in c("ordinal","mice","mitml","ordinalNet","sandwich","lme4","jsonlite"))cat(p,as.character(packageVersion(p)),"\\n")'
    return subprocess.run([str(root/'env/bin/Rscript'),'--vanilla','-e',code],
                           check=True,capture_output=True,text=True).stdout


def code_signature(root: Path, stage: str) -> dict:
    names = {'preflight':['clinical','configuration','reference','cohort'],
             'lesion':['cohort','space','patient','index','spectral','reference'],
             'clinical':['clinical','configuration'], 'annotation':['annotation','spatial_export','plotting']}.get(stage,[])
    paths = [root/'functions'/f'{name}.py' for name in names]
    # 模型适配器决定结局、协变量和输入表；其改变也必须使旧模型失效。
    # lesion 不包含调度模块，因此纯下游修改不会重算病灶特征。
    if stage in ['preflight','clinical','imputation','association','patterns','validation','annotation','report']:
        paths.append(root/'functions/pipeline.py')
    if stage == 'report': paths.append(root/'functions/results_report.py')
    if stage in ['clinical','annotation']: paths.append(root/'stage.py')
    if stage in ['imputation','association','patterns','validation']:
        paths += list((root/'statistics').glob('*.R'))
    signature={str(path.relative_to(root)):digest(path) for path in paths}
    if stage in ['imputation','association','patterns','validation']: signature['r_versions']=r_versions(root)
    if stage in ['clinical','annotation']:
        signature['packages']={name:importlib.metadata.version(name) for name in ['numpy','pandas','pyreadstat','scipy']}
    return signature


def ensure_stage_fresh(root: Path, name: str, config: dict, args, state: dict,
                       verified: set | None = None) -> None:
    """递归核对输入/方法指纹及真实产物；一次检查内避免重复读取同一阶段。"""
    verified = set() if verified is None else verified
    if name in verified: return
    previous = state['stages'].get(name,{})
    if previous.get('status') != 'completed':
        raise ValueError(f'需要先完成阶段 {name}')
    upstream_signatures = {}
    for dependency in dependencies[name]:
        ensure_stage_fresh(root,dependency,config,args,state,verified)
        upstream_signatures[dependency] = state['stages'][dependency]['signature']
    if name == 'association':
        upstream_signatures['clinical_path'] = state['stages']['clinical']['signature']
    expected = fingerprint({'parameters':stage_parameters(name,config),
        'sources':stage_sources(root,name,config,args),'code':code_signature(root,name),
        'upstream':upstream_signatures})
    if expected != previous.get('signature'):
        raise ValueError(f'阶段 {name} 的输入、参数或代码已改变，先续跑该阶段及其下游')
    target = Path(previous['path'])
    complete = read_json(target/'complete.json')
    if complete.get('signature') != expected:
        raise ValueError(f'阶段 {name} 的完成记录与运行指纹不一致')
    verify_files(target,complete['outputs'])
    verified.add(name)


def refresh_report_state(root: Path, config: dict, args, state: dict) -> dict:
    """报告只展示当前且完整的产物；过期或损坏结果保留历史路径并标为 stale。"""
    verified = set()
    for name in stage_order:
        if name == 'report' or state['stages'].get(name,{}).get('status') != 'completed':
            continue
        try:
            ensure_stage_fresh(root,name,config,args,state,verified)
        except (ValueError,KeyError,OSError,subprocess.SubprocessError) as exc:
            state['stages'][name].update(status='stale',reason=f'报告复核未通过：{exc}')
    # 排除心跳、绝对路径及报告自身，避免时间变化或自引用导致指纹漂移。
    return {name:{key:record.get(key) for key in ['status','signature','reason']}
            for name,record in state['stages'].items() if name != 'report'}


def clinical_outcomes(config: dict) -> list[dict]:
    """四个家族先生成完整检验清单；失败模型不从分母删除。"""
    clinical = config['clinical']
    outcomes = []
    for family, definition in config['statistics']['families'].items():
        if not definition.get('enabled',True): continue
        covariates = clinical['covariate_sets'][definition['covariate_set']]
        if family == 'mrs_longitudinal':
            outcomes.append({'name':'mrs_longitudinal','family':family,'type':'longitudinal',
                             'levels':list(range(7)), 'covariates':covariates,
                             'visits':definition['visits'],'visit_columns':list(definition['visits']),
                             'visit_times':list(definition['visits'].values())})
            continue
        for name in definition['outcomes']:
            variable = clinical['variables'][name]
            outcomes.append({'name':name,'family':family,'type':variable['type'],
                             'levels':variable.get('levels'), 'covariates':covariates,
                             'maximum':variable.get('range',[None,None])[-1]})
    return outcomes


def run_child(command: list[str], directory: Path, stop_state: dict, request_path: Path,
              token: str, progress) -> int:
    """只管理本次 Popen 创建的进程组；普通停止交给计算批次检查点。"""
    with (directory/'execution.log').open('a') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        while process.poll() is None:
            if request_path.exists():
                request = read_json(request_path)
                if request.get('token') == token:
                    stop_state['graceful'] = True
                    stop_state['immediate'] |= request.get('immediate',False)
            if stop_state['graceful']:
                (directory/'stop').touch()
                child_control = directory/'control.json'
                if child_control.exists():
                    save_json(directory/'stop_request.json', {'token':read_json(child_control)['token'],
                              'immediate':stop_state['immediate']})
            if stop_state['immediate']:
                # 新建会话的 pgid 等于本子进程 PID；不能终止无关项目。
                try: os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError: pass  # 轮询后子进程恰好退出，无需再次发信号。
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try: os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError: pass
                    process.wait()
                return 130
            progress()
            time.sleep(.5)
        return int(process.returncode)


def statistics_job(root: Path, stage: str, config: dict, upstream: dict, directory: Path) -> dict:
    """R 后端只接受白名单数据与参数，不执行 config 中的表达式。"""
    clinical_dir = Path(upstream.get('clinical') or upstream.get('clinical_path'))
    cohort = 'validation' if stage == 'validation' else 'discovery'
    data = clinical_dir/cohort/'analysis.tsv.gz'
    if not data.exists(): raise FileNotFoundError(f'缺少 {cohort} 临床整合表：{data}')
    outcomes = clinical_outcomes(config)
    if stage == "validation":
        common_file = clinical_dir / "common_variables.json"
        if not common_file.exists():
            raise ValueError("validation 缺少 discovery/validation 共同变量清单，不能绕过共同字段筛选")
        common = set(read_json(common_file).get("variables", []))
        def is_common(outcome):
            required = set(outcome.get("covariates", []))
            if outcome["name"] == "mrs_longitudinal":
                required.update(outcome.get("visit_columns", []))
            else:
                required.add(outcome["name"])
            return required.issubset(common)
        outcomes = [outcome for outcome in outcomes if is_common(outcome)]
        if not outcomes:
            raise ValueError("validation 没有与 discovery 共同的临床结局，跳过该验证阶段")
    return {'stage':stage,'input_csv':str(data), 'output_dir':str(directory),
            'config':config,'variables':config['clinical']['variables'],
            'outcomes':outcomes, 'stop_file':str(directory/'stop'),
            'upstream_imputation_dir':upstream.get('imputation'),
            'frozen_models_dir':upstream.get('patterns'),
            'reference_id':read_json(root/'reference/manifest.json')['reference_id'],
            'feature_columns':[f'{channel}_{mode:03d}' for channel in ['g','w'] for mode in range(1,130)]}


def report_pipeline(directory: Path, state: dict) -> None:
    rows = ''.join(f'<tr><td>{html.escape(name)}</td><td>{html.escape(value.get("status",""))}</td>'
                   f'<td>{html.escape(value.get("reason",""))}</td></tr>' for name,value in state['stages'].items())
    document = '<!doctype html><meta charset="utf-8"><title>正式项目运行报告</title>'
    document += '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto}td,th{padding:10px;border:1px solid #ddd}table{border-collapse:collapse}</style>'
    document += '<h1>LQT135 双通道谐波损伤</h1><p>G：局部组织覆盖的平方载荷损失；W：固定健康尺度的连接支撑损失。脑干为混合组织。</p>'
    document += f'<p>本次状态：{html.escape(state["status"])}</p><table><tr><th>阶段</th><th>状态</th><th>原因</th></tr>{rows}</table>'
    document += '<p>临床输入未配置时不产生临床发现。archive 内的三例与合成结果只用于工程验证。患者形状变化不属于默认流程；不生成患者病灶质控图片。</p>'
    for name,value in state['stages'].items():
        if value.get('path'):
            relative = os.path.relpath(value['path'],directory)
            document += f'<p>{html.escape(name)}：<a href="{relative}/">阶段产物</a></p>'
    from functions.results_report import report_sections
    document += report_sections(directory,state)
    (directory/'report.html').write_text(document,encoding='utf-8')


def run_pipeline(root: Path, args, config: dict) -> int:
    directory = run_directory(root,args.run,config['paths']['results'])
    requested = [part.strip() for part in args.stages.split(',')]
    use_default = requested == ['all']
    if use_default: requested = config['runtime'].get('stages',stage_order[:])
    disabled = set()
    if not config['clinical'].get('enabled',True): disabled.update(stage_order[2:-1])
    if not config['patterns'].get('enabled',True): disabled.update(['patterns','validation','annotation'])
    if not config['validation'].get('enabled',True): disabled.update(['validation','annotation'])
    if not config['annotation'].get('enabled',True): disabled.add('annotation')
    if not use_default and set(requested) & disabled:
        raise ValueError(f'配置禁用了请求阶段：{sorted(set(requested)&disabled)}')
    requested = [name for name in requested if name not in disabled]
    if set(requested)-set(stage_order): raise ValueError('未知阶段名称')
    requested = [name for name in stage_order if name in requested]
    state_path = directory/'stages.json'
    with run_lock(directory):
        if state_path.exists() and not args.resume: raise ValueError('运行已存在，请 --resume 或使用新运行名')
        state = read_json(state_path) if state_path.exists() else {'stages':{}}
        if args.resume and args.subjects is None and state['stages'].get('lesion',{}).get('path'):
            lesion_root = Path(state['stages']['lesion']['path'])
            from functions.cohort import configured_cohorts
            manifests = [lesion_root / cohort / "manifest.json"
                         for cohort in configured_cohorts(config)]
            subject_ids = []
            for prior_manifest in manifests:
                if prior_manifest.exists():
                    subject_ids.extend(item["subject_id"]
                                       for item in read_json(prior_manifest).get("subjects", []))
            if subject_ids:
                args.subjects = ",".join(dict.fromkeys(subject_ids))
        state.update(status='running',heartbeat=time.time())
        save_json(directory/'control.json',identity(root,args.token))
        stopping = {'graceful':False,'immediate':False}
        def on_signal(number,frame):
            stopping['graceful'] = True
            stopping['immediate'] |= number == signal.SIGTERM
        signal.signal(signal.SIGINT,on_signal)
        signal.signal(signal.SIGTERM,on_signal)
        def progress():
            state['heartbeat'] = time.time()
            save_json(state_path,state)
            save_json(directory/'progress.json',state)
        for stage in requested:
            request_path = directory/'stop_request.json'
            if request_path.exists() and read_json(request_path).get('token') == args.token:
                stopping['graceful'] = True
            if stopping['graceful']: break
            record = {'status':'running'}
            try:
                reported_stages = refresh_report_state(root,config,args,state) if stage == 'report' else None
                upstream = {}
                verified = set()
                for dependency in dependencies[stage]:
                    ensure_stage_fresh(root,dependency,config,args,state,verified)
                    previous = state['stages'].get(dependency,{})
                    if previous.get('status') != 'completed':
                        raise ValueError(f'需要先完成阶段 {dependency}')
                    completion = read_json(Path(previous['path'])/'complete.json')
                    verify_files(Path(previous['path']), completion['outputs'])
                    upstream[dependency] = previous['path']
                # association 的临床路径经 MI 清单回溯；不会再次导入或插补。
                if stage == 'association': upstream['clinical_path'] = state['stages']['clinical']['path']
                sources_before=stage_sources(root,stage,config,args)
                code_before=code_signature(root,stage)
                signature_data = {'parameters':stage_parameters(stage,config),
                    'sources':sources_before,'code':code_before,
                    'upstream':{key:read_json(Path(path)/'complete.json')['signature'] for key,path in upstream.items()}}
                if stage == 'report': signature_data['reported_stages'] = reported_stages
                signature = fingerprint(signature_data)
                target = directory/'stages'/stage/signature[:16]
                record.update(path=str(target),signature=signature)
                old_signature = state['stages'].get(stage,{}).get('signature')
                if old_signature and old_signature != signature:
                    invalidated = {stage}
                    for dependent in stage_order:
                        if set(dependencies[dependent]) & invalidated:
                            invalidated.add(dependent)
                            if dependent in state['stages']:
                                state['stages'][dependent]['status'] = 'stale'
                                state['stages'][dependent]['reason'] = '上游结果已改变，需按新指纹重新执行'
                if (target/'complete.json').exists():
                    complete = read_json(target/'complete.json')
                    verify_files(target,complete['outputs'])
                    record['status'] = 'completed'
                    state['stages'][stage] = record
                    progress()
                    continue
                target.mkdir(parents=True,exist_ok=True)
                (target/'stop').unlink(missing_ok=True)
                state['stages'][stage] = record
                progress()
                save_json(target/'parameters.json',stage_parameters(stage,config))
                code = 0
                if stage == 'preflight':
                    from functions.clinical import preflight_clinical
                    from functions.reference import verify_reference,load_reference,numerical_checks
                    from functions.cohort import configured_cohorts
                    verify_reference(root)
                    health = numerical_checks(load_reference(root,with_index=False))
                    issues = {cohort:preflight_clinical(root,config,cohort) for cohort in configured_cohorts(config)}
                    try:
                        stage_sources(root,'lesion',config,args)
                        image_issues = []
                    except (ValueError,FileNotFoundError) as exc: image_issues = [str(exc)]
                    save_json(target/'preflight.json',{'healthy':health,'clinical_issues':issues,
                              'clinical_ready':not any(issues.values()),'imaging_issues':image_issues,
                              'imaging_ready':not image_issues})
                elif stage == 'lesion':
                    # 每个队列使用独立的逐例状态库和结果目录；健康参考和方法指纹相同。
                    cohort_status = {}
                    for cohort in configured_cohorts(config):
                        cohort_target = target / cohort
                        cohort_target.mkdir(parents=True, exist_ok=True)
                        child_config = dict(config, output_directory=str(cohort_target),
                                            active_cohort=cohort)
                        job_config = cohort_target / "job_config.json"
                        save_json(job_config, child_config)
                        command = [sys.executable, "-u", "-B", str(root / "run.py"),
                                   "--config", str(job_config), "--run", args.run,
                                   "--internal-image", "--token", uuid.uuid4().hex]
                        if args.subjects:
                            command += ["--subjects", args.subjects]
                        if (cohort_target / "manifest.json").exists():
                            command += ["--resume"]
                        child_code = run_child(command, cohort_target, stopping,
                                                directory / "stop_request.json",
                                                args.token, progress)
                        cohort_status[cohort] = {"code": child_code}
                        progress_file = cohort_target / "progress.json"
                        if progress_file.exists():
                            cohort_status[cohort]["progress"] = read_json(progress_file)
                        if child_code:
                            code = child_code
                            break
                        if stopping["graceful"]:
                            code = 130
                            break
                    save_json(target / "cohorts.json", cohort_status)
                elif stage in ['clinical','annotation']:
                    if stage=='clinical':
                        from functions.clinical import preflight_clinical
                        required_cohort = config.get("active_cohort", "discovery")
                        issues = preflight_clinical(root, config, required_cohort)
                        if issues:
                            raise ValueError("; ".join(issues))
                    job={'stage':stage,'output_dir':str(target),'config':config,'upstream':upstream,
                         'clinical_dir':state['stages'].get('clinical',{}).get('path')}
                    save_json(target/'job.json',job)
                    code=run_child([sys.executable,'-B',str(root/'stage.py'),str(target/'job.json')],
                                   target,stopping,directory/'stop_request.json',args.token,progress)
                elif stage in ['imputation','association','patterns','validation']:
                    job = statistics_job(root,stage,config,upstream,target)
                    save_json(target/'job.json',job)
                    code = run_child([str(root/'env/bin/Rscript'),'--vanilla',str(root/'statistics/run.R'),
                                      str(target/'job.json')],target,stopping,directory/'stop_request.json',args.token,progress)
                elif stage == 'report': report_pipeline(directory,state)
                backend_status = target/'status.json'
                if backend_status.exists(): record['backend'] = read_json(backend_status)
                if code in [75,130] or stopping['graceful']:
                    record['status'] = 'stopped'
                    stopping['graceful'] = True
                elif code:
                    raise RuntimeError(f'后端退出码 {code}；见 {target}/execution.log')
                else:
                    if sources_before!=stage_sources(root,stage,config,args) or code_before!=code_signature(root,stage):
                        raise ValueError('阶段执行期间输入或代码发生变化，结果不标记完成')
                    excluded = {'complete.json','execution.log','stop','control.json','stop_request.json','progress.json','run.lock','state.sqlite'}
                    outputs = {str(path.relative_to(target)):digest(path) for path in target.rglob('*')
                               if path.is_file() and path.name not in excluded and 'work' not in path.relative_to(target).parts}
                    save_json(target/'complete.json',{'signature':signature,'outputs':outputs})
                    record['status'] = 'completed'
            except Exception as exc:
                record.update(status='blocked' if isinstance(exc,(FileNotFoundError,ValueError)) else 'failed',
                              reason=f'{type(exc).__name__}: {exc}')
                print(f'{stage}: {record["reason"]}',flush=True)
            state['stages'][stage] = record
            progress()
        refresh_report_state(root,config,args,state)
        state['status'] = 'stopped' if stopping['graceful'] else ('completed' if all(
            state['stages'].get(name,{}).get('status') == 'completed' for name in requested) else 'incomplete')
        report_pipeline(directory,state)
        progress()
        return 0 if state['status'] in ['completed','stopped'] else 2
