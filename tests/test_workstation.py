from pathlib import Path
import json
import os
import signal
import subprocess
import sys
import time
import fcntl

import nibabel as nib
import numpy as np
import psutil
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import workstation as ws


@pytest.mark.parametrize('argv,cwd,expected', [
    (['python', '-u', '-B', '/project/run.py', '--token', 'x'], '/elsewhere', 'run'),
    (['python', 'run.py'], '/project', 'run'),
    (['python', 'run.py'], '/other', None),
    (['python', '/other/run.py', '/project/run.py'], '/project', None),
    (['python', '-c', 'print("/project/run.py")'], '/project', None),
    (['python', '/project/monitor.py'], '/other', 'monitor'),
    (['R', '--no-echo', '--file=/project/statistics/run.R'], '/other', 'stage'),
    (['Rscript', '--vanilla', '/project/statistics/run.R', 'job.json'], '/other', 'stage'),
    (['Rscript', '/other/run.R', '/project/statistics/run.R'], '/project', None),
    (['/project/env/bin/python', '-c', 'from multiprocessing.spawn import spawn_main; spawn_main()'], '/project', 'worker'),
    (['/project/env/bin/python', '-c', 'from multiprocessing.spawn import spawn_main; spawn_main()'], '/other', None),
    (['/other/env/bin/python', '-c', 'from multiprocessing.spawn import spawn_main; spawn_main()'], '/project', None),
])
def test_process_scope(argv, cwd, expected):
    assert ws.process_role(Path('/project'), argv, cwd) == expected


def test_signal_does_not_target_recycled_pid(monkeypatch):
    class Process:
        pid = 123
        def create_time(self): return 20
        def terminate(self): raise AssertionError('must not signal reused PID')
    monkeypatch.setattr(ws.psutil, 'Process', lambda pid: Process())
    ws.signal_verified((Process(), 'run', 10))


@pytest.mark.parametrize('ignore_term', [True, False])
def test_stop_tracks_children_and_keeps_other_projects(tmp_path, ignore_term):
    root = tmp_path / 'project'
    other = tmp_path / 'other'
    root.mkdir()
    other.mkdir()
    (root / 'run.py').write_text(
        'import signal, subprocess, sys, time\n'
        f'signal.signal(signal.SIGTERM, signal.{"SIG_IGN" if ignore_term else "SIG_DFL"})\n'
        'child = subprocess.Popen([sys.executable,"-c","import time; time.sleep(60)"])\n'
        'print(child.pid, flush=True)\n'
        'time.sleep(60)\n')
    (other / 'run.py').write_text('import time; time.sleep(60)\n')
    parent = subprocess.Popen([sys.executable, str(root / 'run.py')], cwd=root,
                              stdout=subprocess.PIPE, text=True)
    unrelated = subprocess.Popen([sys.executable, str(other / 'run.py')], cwd=other)
    child_pid = int(parent.stdout.readline().strip())
    try:
        assert unrelated.pid not in ws.inventory(root)
        ws.stop_all(root, wait_seconds=.2)
        parent.wait(timeout=3)
        assert unrelated.poll() is None
        assert not ws.inventory(root)
        if psutil.pid_exists(child_pid):
            assert psutil.Process(child_pid).status() == psutil.STATUS_ZOMBIE
    finally:
        for process in [parent, unrelated]:
            if process.poll() is None:
                process.kill()
            process.wait(timeout=3)
        try:
            psutil.Process(child_pid).kill()
        except psutil.NoSuchProcess:
            pass


def test_orphan_style_spawn_worker_is_visible(tmp_path):
    root = tmp_path / 'project'
    (root / 'env/bin').mkdir(parents=True)
    python = root / 'env/bin/python'
    python.symlink_to(sys.executable)
    worker = subprocess.Popen([str(python), '-c',
        'from multiprocessing.spawn import spawn_main; import time; time.sleep(60)'], cwd=root)
    try:
        assert worker.pid in ws.inventory(root)
        ws.stop_all(root, wait_seconds=.1)
        worker.wait(timeout=3)
    finally:
        if worker.poll() is None:
            worker.kill()
        worker.wait(timeout=3)


def test_held_lock_is_not_removed(tmp_path):
    path = tmp_path / 'results/example/run.lock'
    path.parent.mkdir(parents=True)
    config = {'paths': {'results': 'results'}}
    with path.open('w') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert ws.held_locks(tmp_path, config) == [str(path)]
        with pytest.raises(RuntimeError, match='尚未完全停止'):
            ws.require_stopped(tmp_path, config)
    assert not ws.held_locks(tmp_path, config)
    assert path.exists()


@pytest.fixture
def target(tmp_path):
    root = tmp_path / 'target'
    (root / 'functions').mkdir(parents=True)
    for name in ws.PAYLOAD_FILES:
        (root / name).write_text('# old code\n')
    (root / 'config.workstation.yaml').write_text('private configuration\n')
    (root / 'results').mkdir()
    (root / 'results/private.txt').write_text('private result\n')
    return root


def test_install_verifies_and_backs_up_without_touching_data(target):
    ws.install(target, {'paths': {'results': 'results'}})
    manifest = ws.verified_payload()
    for name, detail in manifest['files'].items():
        assert ws.sha256(target / name) == detail['sha256']
    backups = list((target / 'archive/workstation_updates').iterdir())
    assert len(backups) == 1
    for name in ws.PAYLOAD_FILES:
        assert (backups[0] / name).read_text() == '# old code\n'
    assert (target / 'config.workstation.yaml').read_text() == 'private configuration\n'
    assert (target / 'results/private.txt').read_text() == 'private result\n'


def test_install_refuses_active_processes(target, monkeypatch):
    monkeypatch.setattr(ws, 'inventory', lambda root: {123: (None, 'run', 1)})
    with pytest.raises(RuntimeError, match='尚未完全停止'):
        ws.install(target, {'paths': {'results': 'results'}})
    assert not (target / 'archive').exists()


def test_install_refuses_external_symlink(target, tmp_path):
    outside = tmp_path / 'external.py'
    outside.write_text('# external\n')
    (target / 'run.py').unlink()
    (target / 'run.py').symlink_to(outside)
    with pytest.raises(ValueError, match='位于项目之外'):
        ws.install(target, {'paths': {'results': 'results'}})
    assert outside.read_text() == '# external\n'


def test_install_rolls_back_after_replace_failure(target, monkeypatch):
    replace = ws.os.replace
    calls = []
    def fail_second(source, dest):
        calls.append(dest)
        if len(calls) == 2:
            raise OSError('synthetic disk failure')
        replace(source, dest)
    monkeypatch.setattr(ws.os, 'replace', fail_second)
    with pytest.raises(OSError, match='synthetic disk failure'):
        ws.install(target, {'paths': {'results': 'results'}})
    assert all((target / name).read_text() == '# old code\n' for name in ws.PAYLOAD_FILES)


def test_tampered_manifest_is_rejected(tmp_path):
    (tmp_path / 'manifest.json').write_text(json.dumps({'files': {'../config.yaml': {}}}))
    with pytest.raises(ValueError, match='四个文件'):
        ws.verified_payload(tmp_path)


@pytest.mark.parametrize('scale,shape,valid', [(1, (2, 2, 2), True),
                                             (2, (2, 2, 2), False),
                                             (1, (2, 2, 2, 2), False)])
def test_header_check_is_read_only(tmp_path, scale, shape, valid):
    path = tmp_path / 'mask.nii.gz'
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.uint8), np.diag([scale, scale, scale, 1])), path)
    before = ws.sha256(path)
    if valid:
        assert ws.check_header(path, 1e-4) == 'valid_header'
    else:
        with pytest.raises(ValueError):
            ws.check_header(path, 1e-4)
    assert ws.sha256(path) == before


def test_empty_file_keeps_existing_exclusion_rule(tmp_path):
    path = tmp_path / 'empty.nii.gz'
    path.touch()
    assert ws.check_header(path, 1e-4) == 'empty_file'


def test_maintenance_can_read_config_before_worker_limit_update(tmp_path):
    (tmp_path / 'config.yaml').write_text('paths:\n  results: results\nruntime:\n  workers: 60\n')
    assert ws.maintenance_settings(tmp_path, 'config.yaml')['runtime']['workers'] == 60


@pytest.mark.parametrize('validation_scale', [1, 2])
def test_check_reports_all_headers_and_blocks_invalid_resolution(target, monkeypatch, validation_scale):
    from functions import reference, cohort, clinical
    ws.install(target, {'paths': {'results': 'results'}})
    for name, scale in [('discovery', 1), ('validation', validation_scale)]:
        path = target / f'{name}.nii.gz'
        nib.save(nib.Nifti1Image(np.zeros((2, 2, 2), dtype=np.uint8), np.diag([scale, scale, scale, 1])), path)
    monkeypatch.setattr(reference, 'verify_reference', lambda root: {})
    monkeypatch.setattr(reference, 'load_reference', lambda *args, **kwargs: {})
    monkeypatch.setattr(reference, 'numerical_checks', lambda data: {'synthetic': True})
    monkeypatch.setattr(cohort, 'configured_cohorts', lambda config: ['discovery', 'validation'])
    monkeypatch.setattr(cohort, 'discover', lambda root, selected, config, cohort: [
        {'subject_id': cohort, 'path': str(target / f'{cohort}.nii.gz')}])
    monkeypatch.setattr(clinical, 'preflight_clinical', lambda root, config, name: [])
    config = {'paths': {'results': 'results'}, 'runtime': {'threads': 1, 'stages': []},
              'clinical': {'enabled': True}}
    if validation_scale == 1:
        ws.check(target, config, workers=1)
    else:
        with pytest.raises(RuntimeError, match='检查未通过'):
            ws.check(target, config, workers=1)
    reports = list((target / 'archive/workstation_checks').glob('*.json'))
    assert len(reports) == 1
    report = json.loads(reports[0].read_text())
    assert report['cohorts']['discovery']['total'] == 1
    assert len(report['cohorts']['validation']['header_errors']) == (validation_scale != 1)
