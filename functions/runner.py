"""一、锁定输入；二、有界并行；三、停止和恢复；四、独立计算/人工状态。"""

from __future__ import annotations
import concurrent.futures
import multiprocessing
import os
import signal
import sqlite3
import time
import traceback
from pathlib import Path
import psutil
from functions.storage import read_json, save_json, fingerprint, digest
from functions.reference import verify_reference, load_reference, numerical_checks
from functions.patient import initialize_worker, process_subject, committed_valid
from functions.cohort import discover, method_identity, export_tables
from functions.control import run_directory, run_lock, identity


def run_analysis(root: Path, args, config: dict) -> int:
    directory = Path(config.get("output_directory", run_directory(root, args.run)))
    with run_lock(directory):
        if (directory / "manifest.json").exists() and not args.resume:
            raise ValueError("运行名已存在，续跑请加 --resume，或指定新名称")
        return coordinate(root, directory, args, config)


def coordinate(root: Path, directory: Path, args, config: dict) -> int:
    """唯一状态写入者；进程池仅处理逐例产物。"""
    controls = identity(root, args.token)
    save_json(directory / "control.json", controls)
    stopping = {"graceful": False, "immediate": False}

    def on_signal(number, frame):
        stopping["graceful"] = True
        if number == signal.SIGTERM:
            stopping["immediate"] = True

    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)
    save_json(
        directory / "progress.json",
        {
            "status": "preflight",
            "done": 0,
            "failed": 0,
            "running": 0,
            "queued": 0,
            "heartbeat": time.time(),
        },
    )
    reference = verify_reference(root)
    numerical = numerical_checks(load_reference(root, with_index=False))
    selected = args.subjects
    manifest_path = directory / "manifest.json"
    previous = read_json(manifest_path) if manifest_path.exists() else None
    if args.resume and previous:
        selected = (
            ",".join(item["subject_id"] for item in previous["subjects"])
            if selected is None
            else selected
        )
    selected = selected or config["subjects"].get("default", "all")
    if isinstance(selected, list):
        selected = ",".join(selected)
    cohort = config.get("active_cohort", "discovery")
    subjects = discover(root, selected, config, cohort=cohort)
    method = method_identity(root, config, reference)
    signature = fingerprint(method)
    manifest = {"cohort": cohort, "subjects": subjects, "method": method, "signature": signature}
    if previous:
        if not args.resume:
            raise ValueError("运行名已存在，续跑请加 --resume，或指定新名称")
        if previous != manifest:
            raise ValueError("输入、参考或计算方法指纹改变；禁止混写，请使用新运行名")
    elif args.resume:
        raise ValueError("未找到运行清单，不能续跑")
    else:
        save_json(manifest_path, manifest)
    config = dict(config, reference_id=reference["reference_id"])
    database = sqlite3.connect(directory / "state.sqlite")
    database.execute(
        "CREATE TABLE IF NOT EXISTS jobs (subject TEXT PRIMARY KEY, status TEXT, error TEXT, record_sha256 TEXT)"
    )
    previous_hashes = dict(database.execute("SELECT subject,record_sha256 FROM jobs"))
    for item in subjects:
        target = directory / "subjects" / item["subject_id"]
        status = (
            "done" if committed_valid(target, item["sha256"], signature) else "queued"
        )
        if (
            status == "done"
            and read_json(target / "record.json")["automatic_status"]
            == "excluded_empty_lesion"
        ):
            status = "excluded_empty_lesion"
        if target.exists() and status == "queued":
            raise ValueError(f"已提交结果不完整或被改变：{target}")
        record_hash = digest(target / "record.json") if target.exists() else ""
        if (
            previous_hashes.get(item["subject_id"])
            and previous_hashes[item["subject_id"]] != record_hash
        ):
            raise ValueError(f"已提交病例记录被改变：{target}")
        database.execute(
            "INSERT OR REPLACE INTO jobs VALUES (?,?,?,?)",
            (item["subject_id"], status, "", record_hash),
        )
    database.commit()
    done = {
        row[0]
        for row in database.execute(
            "SELECT subject FROM jobs WHERE status IN ('done','excluded_empty_lesion')"
        )
    }
    pending = [item for item in subjects if item["subject_id"] not in done]
    original_done = len(done)
    active = {}
    started = time.monotonic()
    pool = None
    if pending:
        pool = concurrent.futures.ProcessPoolExecutor(
            max_workers=config["workers"],
            mp_context=multiprocessing.get_context("spawn"),
            initializer=initialize_worker,
            initargs=(str(root), config["threads"], os.getpid()),
        )

    def update_progress(status: str):
        counts = dict(
            database.execute("SELECT status,COUNT(*) FROM jobs GROUP BY status")
        )
        elapsed = time.monotonic() - started
        finished = (
            counts.get("done", 0)
            + counts.get("excluded_empty_lesion", 0)
            - original_done
        )
        rate = finished / elapsed if elapsed and finished else None
        save_json(
            directory / "progress.json",
            {
                "status": status,
                "total": len(subjects),
                "done": counts.get("done", 0),
                "failed": counts.get("failed", 0),
                "running": counts.get("running", 0),
                "queued": counts.get("queued", 0),
                "excluded_empty_lesion": counts.get("excluded_empty_lesion", 0),
                "current_subjects": [item["subject_id"] for item in active.values()],
                "elapsed_seconds": elapsed,
                "subjects_per_minute": rate * 60 if rate else None,
                "eta_seconds": (
                    (len(pending) + len(active)) / rate
                    if rate and finished >= 2
                    else None
                ),
                "heartbeat": time.time(),
                "lesion_qc": "not_performed_user_request",
                "healthy_checks": numerical,
            },
        )

    try:
        while pending or active:
            request = directory / "stop_request.json"
            if request.exists():
                value = read_json(request)
                if value.get("token") == args.token:
                    stopping["graceful"] = True
                    stopping["immediate"] |= value.get("immediate", False)
            if stopping["immediate"]:
                terminated = []
                for process in psutil.Process().children():
                    try:
                        if "spawn_main" in " ".join(process.cmdline()):
                            process.terminate()
                            terminated.append(process)
                    except psutil.Error:
                        pass
                _, alive = psutil.wait_procs(terminated, timeout=3)
                for process in alive:
                    try:
                        process.kill()
                    except psutil.Error:
                        pass
                for item in active.values():
                    database.execute(
                        "UPDATE jobs SET status='queued',error='interrupted' WHERE subject=?",
                        (item["subject_id"],),
                    )
                database.commit()
                break
            while (
                pending and len(active) < config["workers"] and not stopping["graceful"]
            ):
                item = pending.pop(0)
                future = pool.submit(
                    process_subject, str(root), args.run, item, config, signature
                )
                active[future] = item
                database.execute(
                    "UPDATE jobs SET status='running' WHERE subject=?",
                    (item["subject_id"],),
                )
                database.commit()
            update_progress("stopping" if stopping["graceful"] else "running")
            if not active:
                break
            completed, _ = concurrent.futures.wait(
                active, timeout=0.25, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in completed:
                item = active.pop(future)
                try:
                    result = future.result()
                    result_status = (
                        "done"
                        if result["automatic_status"] == "passed"
                        else "excluded_empty_lesion"
                    )
                    record_hash = digest(
                        directory / "subjects" / item["subject_id"] / "record.json"
                    )
                    database.execute(
                        "UPDATE jobs SET status=?,error='',record_sha256=? WHERE subject=?",
                        (result_status, record_hash, item["subject_id"]),
                    )
                    print(
                        f'{item["subject_id"]}: {result_status} ({result["elapsed_seconds"]:.2f}s)',
                        flush=True,
                    )
                except Exception as exc:
                    message = f"{type(exc).__name__}: {exc}"
                    database.execute(
                        "UPDATE jobs SET status='failed',error=? WHERE subject=?",
                        (message, item["subject_id"]),
                    )
                    error_path = directory / "errors" / f'{item["subject_id"]}.txt'
                    error_path.parent.mkdir(exist_ok=True)
                    error_path.write_text(traceback.format_exc(), encoding="utf-8")
                    print(f'{item["subject_id"]}: FAILED {message}', flush=True)
                database.commit()
    finally:
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
    jobs = list(
        database.execute("SELECT subject,status,error FROM jobs ORDER BY subject")
    )
    remaining = any(status in ["queued", "running"] for _, status, _ in jobs)
    failed = any(status == "failed" for _, status, _ in jobs)
    status = (
        "stopped" if remaining else ("completed_with_errors" if failed else "completed")
    )
    if status == "completed" and not any(item[1] == "done" for item in jobs):
        status = "completed_no_analyzable_subjects"
    active.clear()
    update_progress(status)
    export_tables(root, directory, jobs)
    from functions.reporting import write_report

    write_report(root, directory)
    database.close()
    print(f"{args.run}: {status}", flush=True)
    return 2 if failed else 0
