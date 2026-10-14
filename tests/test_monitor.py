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
    assert "discovery" in output and "done=7" in output
    assert "validation" in output and "done=3" in output
    assert "failed=1, running=2, queued=4, excluded_empty_lesion=1" in output

    state["token"] = "previous"
    save_json(tmp_path / "progress.json", state)
    monitor.main()
    output = capsys.readouterr().out
    assert "等待本次运行刷新进度" in output
    assert "done=" not in output
