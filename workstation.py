#!/usr/bin/env python
"""现有项目的停止、补丁安装和启动前检查；不启动患者计算。"""
from __future__ import annotations

import argparse
from datetime import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

import psutil

PACKAGE = Path(__file__).resolve().parent
PAYLOAD_FILES = {'monitor.py', 'run.py', 'functions/pipeline.py',
                 'functions/configuration.py', 'functions/clinical.py'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def project_root(value):
    root = Path(value).resolve(strict=True)
    for name in ['run.py', 'stop.py', 'functions/control.py']:
        if not (root / name).is_file():
            raise ValueError(f'不是现有项目根目录，缺少 {name}: {root}')
    return root


def resolved_argument(value, cwd):
    return (Path(cwd) / value).resolve()


def process_role(root, argv, cwd):
    """按实际入口、或项目环境+cwd+spawn入口识别，不能仅凭命令子串。"""
    if not argv or not cwd:
        return None
    executable = Path(argv[0]).name
    if executable.startswith(('python', 'pypy')):
        index = 1
        while index < len(argv) and argv[index].startswith('-'):
            option = argv[index]
            if option in ['-c', '-m']:
                break
            index += 2 if option in ['-W', '-X'] else 1
        if index < len(argv) and not argv[index].startswith('-'):
            entry = resolved_argument(argv[index], cwd)
            roles = {root / 'run.py': 'run', root / 'stage.py': 'stage',
                     root / 'monitor.py': 'monitor'}
            if entry in roles:
                return roles[entry]
        bootstraps = ('from multiprocessing.spawn import spawn_main',
                      'from multiprocessing.resource_tracker import main')
        # spawn 子进程可能已被收养，保留环境、工作目录和明确启动代码三项约束。
        interpreter = Path(os.path.abspath(os.path.join(cwd, argv[0])))
        if (interpreter.is_relative_to(root / 'env') and Path(cwd).resolve() == root
                and '-c' in argv and any(arg.startswith(bootstraps) for arg in argv)):
            return 'worker'
    if executable in ['R', 'Rscript']:
        entries = [arg.split('=', 1)[1] for arg in argv if arg.startswith('--file=')]
        if executable == 'Rscript':
            entries += [arg for arg in argv[1:] if not arg.startswith('-')][:1]
        if any(resolved_argument(arg, cwd) == root / 'statistics/run.R' for arg in entries):
            return 'stage'
    return None


def inventory(root):
    found = {}
    for proc in psutil.process_iter(['pid', 'uids', 'cmdline', 'cwd', 'status', 'create_time']):
        try:
            info = proc.info
            if (proc.pid == os.getpid() or info['uids'].real != os.getuid()
                    or info['status'] == psutil.STATUS_ZOMBIE):
                continue
            role = process_role(root, info['cmdline'], info['cwd'])
            if role:
                found[proc.pid] = (proc, role, info['create_time'])
        except (psutil.NoSuchProcess, psutil.AccessDenied, OSError):
            continue
    for proc, _, _ in list(found.values()):
        try:
            for child in proc.children(recursive=True):
                if (child.pid != os.getpid() and child.uids().real == os.getuid()
                        and child.status() != psutil.STATUS_ZOMBIE):
                    found.setdefault(child.pid, (child, 'child', child.create_time()))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return found


def print_inventory(found):
    print(f'本用户在本项目的活动进程：{len(found)}', flush=True)
    labels = {'run': '总控/队列', 'stage': '后端阶段', 'monitor': '监控',
              'worker': '工作进程', 'child': '子进程'}
    for pid, (_, role, _) in sorted(found.items()):
        print(f'  PID {pid:8d}  {labels[role]}', flush=True)


def signal_verified(entry, force=False):
    proc, _, created = entry
    try:
        current = psutil.Process(proc.pid)
        if current.create_time() != created or current.uids().real != os.getuid():
            return
        current.kill() if force else current.terminate()
    except psutil.NoSuchProcess:
        pass


def stop_all(root, wait_seconds=30):
    found = inventory(root)
    tracked = dict(found)

    def refresh():
        tracked.update(inventory(root))
        alive = {}
        for pid, entry in tracked.items():
            try:
                proc = psutil.Process(pid)
                if (proc.create_time() == entry[2] and proc.uids().real == os.getuid()
                        and proc.status() != psutil.STATUS_ZOMBIE):
                    alive[pid] = entry
            except psutil.NoSuchProcess:
                continue
        return alive

    print_inventory(found)
    deadline = time.monotonic() + wait_seconds
    sent = set()
    while found and time.monotonic() < deadline:
        for pid, entry in found.items():
            identity = (pid, entry[2])
            # 先让协调进程执行自身的退出清理；孤立 worker 在超时后单独处理。
            if entry[1] not in ['worker', 'child'] and identity not in sent:
                signal_verified(entry)
                sent.add(identity)
        time.sleep(.25)
        found = refresh()
    if found:
        print(f'退出等待结束，终止仍存活的 {len(found)} 个已核验项目进程。', flush=True)
        for entry in reversed(list(found.values())):
            signal_verified(entry, force=True)
        deadline = time.monotonic() + 5
        while found and time.monotonic() < deadline:
            time.sleep(.25)
            found = refresh()
    print_inventory(found)
    if found:
        raise RuntimeError('仍有项目进程，停止未完成；不要安装或重启。')


def settings(root, config_file):
    sys.path.insert(0, str(root))
    from functions.configuration import load_config
    return load_config(root / config_file)


def maintenance_settings(root, config_file):
    # 旧版本可能尚不接受 workers=60；停止和安装只需读取结果路径。
    import yaml
    config = yaml.safe_load((root / config_file).read_text(encoding='utf-8'))
    if (not isinstance(config, dict) or not isinstance(config.get('paths'), dict)
            or not isinstance(config['paths'].get('results'), str)
            or not config['paths']['results']):
        raise ValueError('配置必须提供非空 paths.results，才能核查运行锁')
    return config


def held_locks(root, config):
    base = (root / config['paths']['results']).resolve()
    locks = list(base.glob('*/run.lock')) + list(base.glob('*/stages/lesion/*/*/run.lock'))
    held = []
    for path in locks:
        try:
            with path.open('r+') as handle:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    held.append(str(path))
                else:
                    fcntl.flock(handle, fcntl.LOCK_UN)
        except FileNotFoundError:
            continue
    return held


def require_stopped(root, config):
    found = inventory(root)
    print_inventory(found)
    locks = held_locks(root, config)
    if found or locks:
        raise RuntimeError(f'项目尚未完全停止；仍持有的运行锁：{locks}')
    print('停止检查通过：活动进程 0，持有的运行锁 0。', flush=True)


def verified_payload(package=PACKAGE):
    manifest = json.loads((package / 'manifest.json').read_text())
    if set(manifest['files']) != PAYLOAD_FILES:
        raise ValueError('补丁文件列表与允许替换的代码文件不符')
    for name, detail in manifest['files'].items():
        path = package / 'payload' / name
        if not path.resolve().is_relative_to((package / 'payload').resolve()):
            raise ValueError(f'补丁文件路径越界：{name}')
        if sha256(path) != detail['sha256']:
            raise ValueError(f'补丁校验失败：{name}')
        compile(path.read_bytes(), str(path), 'exec')
    return manifest


def install(root, config, package=PACKAGE):
    manifest = verified_payload(package)
    require_stopped(root, config)
    for name in manifest['files']:
        target = root / name
        if not target.is_file() or not target.resolve().is_relative_to(root):
            raise ValueError(f'安装目标不存在或位于项目之外：{target}')
    base = root / 'archive/workstation_updates'
    if not base.resolve().is_relative_to(root):
        raise ValueError('代码备份目录不能位于项目之外')
    base.mkdir(parents=True, exist_ok=True)
    backup = Path(tempfile.mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'), dir=base))
    for name in manifest['files']:
        saved = backup / name
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / name, saved)
    # 所有原文件备份完成后再写；单文件替换为原子操作，异常时恢复已替换文件。
    replaced = []
    try:
        for name in manifest['files']:
            target = root / name
            descriptor, temporary = tempfile.mkstemp(prefix=f'.{target.name}.update-', dir=target.parent)
            os.close(descriptor)
            try:
                shutil.copy2(package / 'payload' / name, temporary)
                os.replace(temporary, target)
                replaced.append(name)
            finally:
                Path(temporary).unlink(missing_ok=True)
        for name, detail in manifest['files'].items():
            if sha256(root / name) != detail['sha256']:
                raise RuntimeError(f'安装后校验失败：{name}')
    except Exception:
        for name in replaced:
            shutil.copy2(backup / name, root / name)
        raise
    (backup / 'receipt.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'更新完成，旧代码备份：{backup}', flush=True)
    print('配置、患者数据、参考资源、环境和已有结果均保留。', flush=True)


def check_header(path, tolerance):
    import nibabel as nib
    import numpy as np
    if path.stat().st_size == 0:
        return 'empty_file'
    image = nib.load(path)
    if len(image.shape) != 3:
        raise ValueError('病灶不是三维图像')
    if (not np.isfinite(image.affine).all()
            or abs(np.linalg.det(image.affine[:3, :3])) < 1e-12):
        raise ValueError('病灶仿射矩阵无效')
    sizes = nib.affines.voxel_sizes(image.affine)
    if not np.allclose(sizes, 1.0, atol=tolerance, rtol=0):
        raise ValueError(f'病灶体素必须为 1 mm，实际为 {sizes.tolist()}')
    return 'valid_header'


def check(root, config, workers, package=PACKAGE):
    require_stopped(root, config)
    manifest = verified_payload(package)
    for name, detail in manifest['files'].items():
        if sha256(root / name) != detail['sha256']:
            raise ValueError(f'代码与本补丁不一致，请先 install：{name}')
    if not 1 <= workers <= 64:
        raise ValueError('workers 必须在 1—64 之间')
    available = len(os.sched_getaffinity(0))
    if workers > available:
        raise ValueError(f'请求 {workers} 个进程，但当前 CPU 亲和范围只有 {available} 个 CPU')
    if int(config['runtime']['threads']) != 1:
        raise ValueError('本次多进程方案要求配置 runtime.threads: 1')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
        os.environ[name] = '1'
    from functions.reference import verify_reference, load_reference, numerical_checks
    from functions.cohort import discover, configured_cohorts
    from functions.clinical import preflight_clinical, assert_disjoint
    from functions.pipeline import r_versions, stage_order
    report = {'source_commit': manifest['source_commit'], 'workers': workers,
              'available_cpus': available, 'available_memory_gib': round(psutil.virtual_memory().available / 2**30, 1),
              'issues': [], 'cohorts': {}, 'scope': 'reference, config, IDs, SAS fields and NIfTI headers'}
    print('检查参考及数值定义……', flush=True)
    verify_reference(root)
    report['healthy'] = numerical_checks(load_reference(root, with_index=False))
    cohorts = {}
    for name in configured_cohorts(config):
        print(f'{name}: 检查文件清单和 ID（需要读取病灶文件计算校验值）……', flush=True)
        items = discover(root, 'all', config, cohort=name)
        cohorts[name] = {item['subject_id'] for item in items}
        errors, empty = [], 0
        tolerance = config.get('imaging', {}).get('voxel_tolerance_mm', 1e-4)
        for index, item in enumerate(items, 1):
            try:
                empty += check_header(root / item['path'], tolerance) == 'empty_file'
            except Exception as error:
                errors.append({'subject_id': item['subject_id'], 'error': str(error)})
            if index % 1000 == 0 or index == len(items):
                print(f'{name}: 图像头检查 {index}/{len(items)}，异常 {len(errors)}', flush=True)
        clinical = preflight_clinical(root, config, name) if config['clinical'].get('enabled', True) else []
        report['cohorts'][name] = {'total': len(items), 'zero_byte_files': empty,
                                   'header_errors': errors, 'clinical_issues': clinical}
        report['issues'].extend(clinical)
        if errors:
            report['issues'].append(f'{name}: {len(errors)} 个图像头异常')
    if 'validation' in cohorts:
        assert_disjoint(cohorts['discovery'], cohorts['validation'])
    stages = config['runtime'].get('stages', stage_order)
    if config['clinical'].get('enabled', True) and set(stages) & {'imputation', 'association', 'patterns', 'validation'}:
        report['r_versions'] = r_versions(root)
    require_stopped(root, config)
    output_dir = root / 'archive/workstation_checks'
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / (datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '.json')
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'检查记录：{output}', flush=True)
    print('图像头检查不评估解剖配准；二值性由正式病灶阶段检查，配准须由输入处理记录确认。', flush=True)
    if report['issues']:
        for issue in report['issues'][:20]:
            print(issue, flush=True)
        raise RuntimeError('检查未通过，请按检查记录修正后重试。')
    print('CHECK PASSED：可以使用原运行名 --resume 启动所有已启用阶段。', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['status', 'stop', 'install', 'check'])
    parser.add_argument('--project', required=True)
    parser.add_argument('--config', default='config.workstation.yaml')
    parser.add_argument('--workers', type=int, default=60)
    parser.add_argument('--wait-seconds', type=int, default=30)
    args = parser.parse_args()
    if not 1 <= args.wait_seconds <= 300:
        parser.error('--wait-seconds 必须为 1—300')
    root = project_root(args.project)
    if args.action == 'stop':
        stop_all(root, args.wait_seconds)
    config = (settings(root, args.config) if args.action == 'check'
              else maintenance_settings(root, args.config))
    if args.action in ['status', 'stop']:
        require_stopped(root, config)
    elif args.action == 'install':
        install(root, config)
    else:
        check(root, config, args.workers)


if __name__ == '__main__':
    try:
        main()
    except (Exception, KeyboardInterrupt) as error:
        print(f'未完成：{type(error).__name__}: {error}', file=sys.stderr, flush=True)
        raise SystemExit(2)
