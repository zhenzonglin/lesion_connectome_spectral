"""统计后端的 JSON 入口、恢复与假设诊断回归测试（仅合成数据）。"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def statistical_job() -> tuple[Path, Path, dict]:
    """合成 180 例，测试产物只能进入项目 archive。"""
    root = Path(__file__).resolve().parents[1]
    executable = root / "env/bin/Rscript"
    if not executable.exists():
        pytest.skip("项目 R 环境尚未安装")
    parent = root / "archive/tests"
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="statistics_", dir=parent))
    generator = np.random.default_rng(5291)
    size = 180
    frame = pd.DataFrame({
        "subject_id": [f"{index:04d}" for index in range(size)],
        "age": generator.normal(60, 10, size),
        "sex": generator.integers(1, 3, size),
        "log_lesion_volume": generator.normal(8, 1, size),
    })
    frame.loc[0, "subject_id"] = "NA"
    features = [f"{channel}_{mode:03d}" for channel in ("g", "w") for mode in (1, 2)]
    for feature in features:
        frame[feature] = generator.uniform(0, 1, size)
    latent = 3 * frame.g_001 + 2 * frame.w_001 + generator.logistic(size=size)
    frame["mrs_3m"] = np.digitize(latent, [1, 3])
    frame["moca_3m"] = 29 - 3 * frame.g_001 - 2 * frame.w_001 + generator.normal(size=size)
    input_path = directory / "synthetic.tsv.gz"
    frame.to_csv(input_path, sep="\t", index=False)
    covariates = ["age", "sex", "log_lesion_volume"]
    job = {
        "stage": "association", "input_csv": str(input_path),
        "output_dir": str(directory / "association"), "feature_columns": features,
        "mode_pairs": [{"mode": mode, "g": f"g_{mode:03d}", "w": f"w_{mode:03d}"} for mode in (1, 2)],
        "variables": {"age": {"type": "continuous"}, "sex": {"type": "binary", "levels": [1, 2]},
                      "log_lesion_volume": {"type": "continuous"}},
        "outcomes": [
            {"name": "mrs_3m", "type": "ordinal", "family": "mrs_single", "levels": [0, 1, 2], "covariates": covariates},
            {"name": "moca_3m", "type": "continuous", "family": "moca", "maximum": 30, "covariates": covariates},
        ],
        "config": {"statistics": {"mi": {"m": 3, "maxit": 4, "increment": 2, "max_iterations": 4}}},
    }
    return root, directory, job


def invoke(root: Path, directory: Path, job: dict, name: str = "job") -> subprocess.CompletedProcess:
    """通过公开 JSON 接口调用 R，不直接加载内部函数。"""
    job_path = directory / f"{name}.json"
    job_path.write_text(json.dumps(job), encoding="utf-8")
    return subprocess.run([str(root / "env/bin/Rscript"), str(root / "statistics/run.R"), str(job_path)],
                          capture_output=True, text=True, timeout=180, check=False)


def test_r_json_entry_resume_and_parameter_guard(statistical_job: tuple) -> None:
    """实际 CLI、压缩 TSV、ID、比例优势诊断、恢复及参数变化一起核验。"""
    root, directory, job = statistical_job
    first = invoke(root, directory, job)
    assert first.returncode == 0, first.stderr
    output = Path(job["output_dir"])
    table = pd.read_csv(output / "associations.csv")
    assert len(table) == 4 and table.p_joint.notna().all()
    diagnostics = pd.read_csv(output / "model_diagnostics.csv")
    ordinal = diagnostics[diagnostics.outcome == "mrs_3m"]
    assert ordinal.po_min_p.notna().all(), "动态公式或局部数据不得导致比例优势检验失效"
    participants = pd.read_csv(output / "imputation/mrs_3m/participants.csv", dtype=str, keep_default_na=False)
    assert participants.subject_id.iloc[:2].tolist() == ["NA", "0001"]
    checkpoint = output / "mrs_3m/mode_001.rds"
    before = (checkpoint.read_bytes(), checkpoint.stat().st_mtime_ns)
    second = invoke(root, directory, job)
    assert second.returncode == 0, second.stderr
    assert before == (checkpoint.read_bytes(), checkpoint.stat().st_mtime_ns)
    changed = json.loads(json.dumps(job))
    changed["config"]["statistics"]["fdr"] = {"q": 0.1}
    rejected = invoke(root, directory, changed, "changed")
    assert rejected.returncode != 0 and "已改变" in rejected.stderr


def test_r_safe_stop_boundary(statistical_job: tuple) -> None:
    """已请求停止时不派发模型，并保留清楚的 stopped 状态。"""
    root, directory, original = statistical_job
    job = json.loads(json.dumps(original))
    job["output_dir"] = str(directory / "stopped")
    marker = directory / "stop_request"
    marker.write_text("stop", encoding="utf-8")
    job["stop_file"] = str(marker)
    result = invoke(root, directory, job, "stop")
    assert result.returncode == 75
    status = json.loads((Path(job["output_dir"]) / "status.json").read_text(encoding="utf-8"))
    assert status["status"] == "stopped"
    assert not (Path(job["output_dir"]) / "associations.csv").exists()
