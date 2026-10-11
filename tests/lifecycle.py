"""工程验收脚本：在项目自有三例副本上测试停止、锁、指纹与提交恢复。"""

from pathlib import Path
import json
import os
import subprocess
import sys
import time
import uuid
import hashlib
import shutil
import nibabel as nib
import numpy as np
import psutil

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from functions.storage import read_json, digest, save_json
from functions.control import request_stop, matching_process


def wait_for(predicate, limit=180):
    started = time.monotonic()
    while time.monotonic() - started < limit:
        if predicate():
            return
        time.sleep(0.05)
    raise TimeoutError("生命周期测试等待超时")


def launch(name, extra=None):
    command = [
        sys.executable,
        "-u",
        "-B",
        str(root / "run.py"),
        "--run",
        name,
        "--workers",
        "1",
        *(extra or []),
    ]
    log = (root / "checks" / f"{name}.log").open("a")
    process = subprocess.Popen(
        command, cwd=root, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
    )
    log.close()
    return process


def state(name):
    path = root / "results" / name / "progress.json"
    return read_json(path) if path.exists() else {}


def hashes(name):
    return {
        str(p.relative_to(root)): digest(p)
        for p in (root / "results" / name / "subjects").rglob("*")
        if p.is_file()
    }


def main():
    report = {}
    first = "check_graceful"
    process = launch(first)
    wait_for(lambda: state(first).get("running", 0) > 0)
    parent = psutil.Process(read_json(root / "results" / first / "control.json")["pid"])
    wait_for(
        lambda: any(
            "spawn_main" in " ".join(child.cmdline()) for child in parent.children()
        )
    )
    held = [
        child
        for child in parent.children()
        if "spawn_main" in " ".join(child.cmdline())
    ]
    for child in held:
        child.suspend()
    duplicate = subprocess.run(
        [sys.executable, str(root / "run.py"), "--run", first, "--resume"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    assert duplicate.returncode != 0 and "重复启动" in duplicate.stderr
    report["duplicate_lock"] = True
    request_stop(root, first)
    wait_for(lambda: state(first).get("status") == "stopping")
    for child in held:
        child.resume()
    assert process.wait(timeout=180) == 0
    assert state(first)["status"] == "stopped" and state(first)["done"] >= 1
    before = hashes(first)
    manual = root / "qc/user_notes.txt"
    assert not manual.exists()
    manual.write_text("用户方法审核备注，自动运行不得修改。")
    manual_hash = digest(manual)
    process = launch(first, ["--resume"])
    assert process.wait(timeout=180) == 0 and state(first)["status"] == "completed"
    assert all(digest(root / path) == value for path, value in before.items())
    assert digest(manual) == manual_hash
    report["graceful_resume_completed_unchanged"] = True
    report["manual_review_preserved"] = True
    manual.unlink()
    second = "check_immediate"
    process = launch(second)
    work = root / "results" / second / "work"
    wait_for(lambda: state(second).get("running", 0) > 0)
    request_stop(root, second, True)
    assert process.wait(timeout=180) == 0 and state(second)["status"] == "stopped"
    before = hashes(second)
    # 明确模拟一次未提交的截断文件，不伪装成真实患者结果。
    partial = work / "TMS001" / "result.npz"
    if not (root / "results" / second / "subjects/TMS001").exists():
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.write_bytes(b"incomplete_archive")
    process = launch(second, ["--resume"])
    assert process.wait(timeout=180) == 0 and state(second)["status"] == "completed"
    assert all(digest(root / path) == value for path, value in before.items())
    assert not list(work.glob("*/result.npz"))
    report["immediate_stop_and_simulated_partial_write_resume"] = True
    fake = read_json(root / "results" / second / "control.json")
    fake["pid"] = os.getpid()
    fake["created"] = 0
    assert matching_process(root, fake) is None
    report["stale_pid_rejected"] = True
    source = next((root / "lesion").glob("sub-TMS001*.nii.gz"))
    test = (
        root
        / "lesion/sub-TESTCASEses01_space-MNI152NLin6Asym_res-01_label-lesion_mask.nii.gz"
    )
    assert not test.exists()
    shutil.copy2(source, test)
    try:
        third = "check_input"
        process = launch(third, ["--subjects", "TESTCASE"])
        assert process.wait(timeout=180) == 0
        unchanged = hashes(third)
        image = nib.load(test)
        mask = np.asarray(image.dataobj).copy()
        mask[90, 90, 90] = 1 - mask[90, 90, 90]
        nib.save(nib.Nifti1Image(mask, image.affine), test)
        process = launch(third, ["--resume"])
        assert process.wait(timeout=180) != 0
        assert unchanged == hashes(third)
        report["changed_input_blocks_resume"] = True
        nib.save(
            nib.Nifti1Image(np.zeros(image.shape, dtype=np.uint8), image.affine), test
        )
        process = launch("check_empty", ["--subjects", "TESTCASE"])
        assert process.wait(timeout=180) == 0
        assert state("check_empty")["status"] == "completed_no_analyzable_subjects"
        assert state("check_empty")["excluded_empty_lesion"] == 1
        empty_directory = root / "results/check_empty/subjects/TESTCASE"
        assert {p.name for p in empty_directory.iterdir()} == {"record.json"}
        before_empty = digest(empty_directory / "record.json")
        process = launch("check_empty", ["--resume"])
        assert (
            process.wait(timeout=180) == 0
            and digest(empty_directory / "record.json") == before_empty
        )
        import pandas as pd

        assert (
            len(pd.read_csv(root / "results/check_empty/excluded.tsv", sep="\t")) == 1
        )
        assert (
            len(
                pd.read_csv(root / "results/check_empty/mode_features.tsv.gz", sep="\t")
            )
            == 0
        )
        report["empty_recorded_excluded_not_zero_and_resume_preserved"] = True
        bad = np.zeros(image.shape, dtype=np.float32)
        bad[90, 90, 90] = 0.5
        nib.save(nib.Nifti1Image(bad, image.affine), test)
        process = launch("check_failure", ["--subjects", "TESTCASE"])
        assert process.wait(timeout=180) != 0
        assert state("check_failure")["status"] == "completed_with_errors"
        assert not (
            root / "results/check_failure/subjects/TESTCASE/result.npz"
        ).exists()
        report["failure_never_becomes_zero_result"] = True
    finally:
        test.unlink()
    report["passed"] = True
    assert not list((root / "results").glob("*/subjects/*/qc.png"))
    assert not (root / "qc/manual_review.tsv").exists()
    report["no_patient_lesion_qc_generated"] = True
    save_json(root / "checks/lifecycle.json", report)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
