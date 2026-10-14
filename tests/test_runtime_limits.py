"""资源配置通过真实 CLI 校验；不启动患者计算或大量工作进程。"""
import os
import sys
from pathlib import Path

import pytest
import yaml

import run as entrypoint
from functions import pipeline


@pytest.mark.parametrize("configured,override,expected", [
    (60, None, 60),
    (4, 60, 60),
    (65, None, None),
    (4, 65, None),
])
def test_worker_limit_from_config_and_cli(tmp_path, monkeypatch, configured, override, expected):
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))
    config["runtime"].update(workers=configured, threads=1)
    path = tmp_path / "runtime.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    argv = ["run.py", "--config", str(path), "--run", "runtime_test", "--token", "test-token"]
    if override is not None:
        argv += ["--workers", str(override)]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(os, "environ", os.environ.copy())
    dispatched = []

    def run_pipeline(root, args, config):
        dispatched.append(config["workers"])
        assert config["runtime"]["workers"] == expected
        assert config["threads"] == 1
        assert os.environ["OMP_NUM_THREADS"] == "1"
        return 0

    monkeypatch.setattr(pipeline, "run_pipeline", run_pipeline)
    if expected is None:
        with pytest.raises(ValueError, match="64"):
            entrypoint.main()
        assert not dispatched
    else:
        assert entrypoint.main() == 0
        assert dispatched == [expected]
