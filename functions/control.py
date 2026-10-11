"""运行锁、进程身份核查和停止请求；不使用全局 kill 或 WSL 关闭操作。"""

from __future__ import annotations
import contextlib
import fcntl
import os
import re
import time
from pathlib import Path
import psutil
from functions.storage import read_json, save_json


def run_directory(root: Path, name: str, base: str = "results") -> Path:
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", name):
        raise ValueError("运行名仅允许英文、数字、下划线或短横线")
    return (root / base).resolve() / name


@contextlib.contextmanager
def run_lock(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "run.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("相同运行正在执行，拒绝重复启动") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def identity(root: Path, token: str) -> dict:
    process = psutil.Process()
    return {
        "pid": process.pid,
        "created": process.create_time(),
        "root": str(root.resolve()),
        "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
        "token": token,
        "started": time.time(),
    }


def matching_process(root: Path, record: dict):
    """PID、创建时间、启动标识、项目路径和系统启动 ID 均匹配才允许发信号。"""
    try:
        if record["root"] != str(root.resolve()):
            return None
        if (
            record["boot_id"]
            != Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        ):
            return None
        process = psutil.Process(record["pid"])
        command = process.cmdline()
        if abs(process.create_time() - record["created"]) > 0.01:
            return None
        if record["token"] not in command or str(root / "run.py") not in command:
            return None
        if process.status() == psutil.STATUS_ZOMBIE:
            return None
        return process
    except (KeyError, OSError, psutil.Error):
        return None


def request_stop(root: Path, name: str, immediate: bool = False, base: str = "results") -> str:
    directory = run_directory(root, name, base)
    record = read_json(directory / "control.json")
    process = matching_process(root, record)
    if process is None:
        return "未发现身份匹配的活动任务；没有向任何进程发送信号。"
    save_json(
        directory / "stop_request.json",
        {"token": record["token"], "immediate": immediate, "requested": time.time()},
    )
    # 常规停止只写请求；立即停止也交给协调进程处理其自有工作进程。
    if immediate:
        process.terminate()
    return (
        "已请求立即停止，未提交病例将在续跑时重算。"
        if immediate
        else "已停止新病例派发，正在执行的病例提交后退出。"
    )
