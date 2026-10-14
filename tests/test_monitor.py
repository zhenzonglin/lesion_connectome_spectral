"""续跑监控不能把前一次终止状态当成本次运行失败。"""
from copy import deepcopy
import sys

import pytest

import monitor
from functions import configuration
from functions.storage import save_json


@pytest.mark.parametrize("arguments,watch,interval", [
    ([], True, 30),
    (["--watch"], True, 30),
    (["--once"], False, 30),
    (["--interval", "5"], True, 5),
])
def test_monitor_handles_previous_progress_while_new_process_is_alive(
    tmp_path, monkeypatch, capsys, arguments, watch, interval
):
    snapshots = [
        {"status": "incomplete", "heartbeat": 10,
         "stages": {"lesion": {"status": "failed", "reason": "OLD_ERROR"}}},
        {"status": "completed", "heartbeat": 30, "stages": {}},
    ]
    state = {"index": 0}
    monkeypatch.setattr(sys, "argv", ["monitor.py", "--run", "test"] + arguments)
    monkeypatch.setattr(configuration, "load_config", lambda path: {
        "runtime": {"monitor_interval_seconds": 30}, "paths": {"results": "results"},
    })
    monkeypatch.setattr(monitor, "run_directory", lambda *args: tmp_path)
    monkeypatch.setattr(monitor, "read_json", lambda path: (
        deepcopy(snapshots[state["index"]]) if path.name == "progress.json"
        else {"pid": 123, "started": 20}
    ))
    monkeypatch.setattr(monitor, "matching_process", lambda *args: object() if state["index"] == 0 else None)

    sleeps = []
    def advance(seconds):
        sleeps.append(seconds)
        if state["index"] == 0:
            state["index"] += 1
        else:
            raise KeyboardInterrupt

    monkeypatch.setattr(monitor.time, "sleep", advance)
    if watch:
        with pytest.raises(KeyboardInterrupt):
            monitor.main()
    else:
        monitor.main()
    output = capsys.readouterr().out
    assert "OLD_ERROR" not in output
    assert "等待本次运行刷新进度" in output
    assert state["index"] == int(watch)
    assert sleeps == ([interval, interval] if watch else [])
    if watch:
        assert "阶段进度 | completed | 进程存活 False" in output
        assert "继续等待状态更新" in output


@pytest.mark.parametrize("interval", ["0", "-1", "0.5", "nan", "inf"])
def test_monitor_rejects_invalid_intervals(monkeypatch, interval):
    monkeypatch.setattr(sys, "argv", ["monitor.py", "--run", "test", "--interval", interval])
    monkeypatch.setattr(configuration, "load_config", lambda path: {
        "runtime": {"monitor_interval_seconds": 30}, "paths": {"results": "results"},
    })
    with pytest.raises(ValueError, match="有限数值"):
        monitor.main()


def test_monitor_waits_for_missing_run_without_creating_files(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["monitor.py", "--run", "test"])
    monkeypatch.setattr(configuration, "load_config", lambda path: {
        "runtime": {"monitor_interval_seconds": 30}, "paths": {"results": "results"},
    })
    monkeypatch.setattr(monitor, "run_directory", lambda *args: tmp_path)
    def interrupt(seconds):
        assert seconds == 30
        raise KeyboardInterrupt
    monkeypatch.setattr(monitor.time, "sleep", interrupt)
    with pytest.raises(KeyboardInterrupt):
        monitor.main()
    assert "等待总控身份记录" in capsys.readouterr().out
    assert not list(tmp_path.iterdir())


def test_monitor_displays_both_cohorts_and_rejects_old_token(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["monitor.py", "--run", "test", "--once"])
    monkeypatch.setattr(configuration, "load_config", lambda path: {
        "runtime": {"monitor_interval_seconds": 30}, "paths": {"results": "results"},
    })
    monkeypatch.setattr(monitor, "run_directory", lambda *args: tmp_path)
    monkeypatch.setattr(monitor, "matching_process", lambda *args: object())
    stage_path = tmp_path / "stages/lesion/current"
    for cohort, done in [("discovery", 7), ("validation", 3)]:
        save_json(stage_path / cohort / "progress.json", {
            "status": "running", "done": done, "failed": 1, "running": 2,
            "queued": 4, "excluded_empty_lesion": 1,
        })
    state = {"status": "running", "heartbeat": 30, "token": "current",
             "stages": {"lesion": {"status": "running", "path": str(stage_path)}}}
    save_json(tmp_path / "progress.json", state)
    save_json(tmp_path / "control.json", {"pid": 123, "started": 20, "token": "current"})
    monitor.main()
    output = capsys.readouterr().out
    rows = {parts[0]: parts for line in output.splitlines() if ' | ' in line
            for parts in [[part.strip() for part in line.split('|')]]}
    assert rows['阶段/队列'][2:] == ['已完成', '失败', '进行中', '待完成', '排除', '总人数']
    assert rows['discovery'][2:] == ['7', '1', '2', '4', '1', '15']
    assert rows['validation'][2:] == ['3', '1', '2', '4', '1', '11']
    assert rows['lesion / 合计'][2:] == ['10', '2', '4', '8', '2', '26']

    state["token"] = "previous"
    save_json(tmp_path / "progress.json", state)
    monitor.main()
    output = capsys.readouterr().out
    assert "等待本次运行刷新进度" in output
    assert "阶段/队列" not in output


def test_planned_cohort_is_queued_before_it_starts(tmp_path):
    path = tmp_path / 'lesion'
    stage = {'status': 'running', 'path': str(path),
             'cohort_totals': {'discovery': 5, 'validation': 7}}
    save_json(path / 'discovery/progress.json', {
        'status': 'completed', 'total': 5, 'done': 3, 'failed': 1,
        'running': 0, 'queued': 0, 'excluded_empty_lesion': 1,
    })
    rows = monitor.lesion_rows(stage, ['discovery', 'validation'])
    assert rows[2][2] == {'total': 7, 'done': 0, 'failed': 0, 'running': 0,
                          'queued': 7, 'excluded_empty_lesion': 0}
    total = rows[0][2]
    assert total['total'] == sum(total[key] for key in monitor.COUNT_KEYS[:-1]) == 12
    assert total['done'] == 3 and total['queued'] == 7


def test_unknown_or_revalidating_counts_are_not_zero(tmp_path):
    stage = {'status': 'running', 'path': str(tmp_path),
             'cohort_totals': {'discovery': 5, 'validation': 7}}
    save_json(tmp_path / 'discovery/manifest.json', {'subjects': []})
    save_json(tmp_path / 'discovery/progress.json', {
        'status': 'preflight', 'done': 0, 'failed': 0, 'running': 0, 'queued': 0,
    })
    rows = monitor.lesion_rows(stage, ['discovery', 'validation'])
    assert rows[1][2]['done'] is None
    assert rows[1][2]['total'] == 5
    assert rows[0][2]['done'] is None and rows[0][2]['total'] == 12
    assert monitor.patient_counts({}) == dict.fromkeys(monitor.COUNT_KEYS)


def test_legacy_counts_and_nonpatient_stage_use_correct_units(tmp_path, capsys):
    save_json(tmp_path / 'progress.json', {
        'status': 'completed_with_errors', 'done': 4, 'failed': 2,
        'running': 0, 'queued': 0, 'excluded_empty_lesion': 1,
    })
    stages = {
        'lesion': {'status': 'failed', 'path': str(tmp_path)},
        'association': {'status': 'blocked', 'reason': 'synthetic model diagnostic'},
    }
    monitor.print_stage_table(stages, ['discovery', 'validation'])
    output = capsys.readouterr().out
    rows = {parts[0]: parts for line in output.splitlines() if ' | ' in line
            for parts in [[part.strip() for part in line.split('|')]]}
    assert rows['lesion / 合计'][2:] == ['4', '2', '0', '0', '1', '7']
    assert rows['association'][2:] == ['—'] * 6
    assert 'discovery' not in rows and 'validation' not in rows
    assert 'synthetic model diagnostic' in output
