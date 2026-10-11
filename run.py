#!/usr/bin/env python
"""正式总控：启动、限定被试、后台运行和校验后续跑。"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path


def main() -> int:
    # 一、入口参数；环境变量必须在数值库导入前设定。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--run", required=True)
    parser.add_argument(
        "--subjects", help="逗号分隔患者号；all 为全部；续跑省略则使用原清单"
    )
    parser.add_argument("--cohort", choices=["discovery", "validation"],
                        help="直接运行单个队列；正式总控默认同时运行双队列")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--stages", default="all", help="逗号分隔阶段；默认 all，单独病灶使用 preflight,lesion")
    parser.add_argument("--preflight", action="store_true", help="只核查配置与参考，不启动患者计算")
    parser.add_argument("--internal-image", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--detach", action="store_true")
    parser.add_argument("--token", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    import yaml
    from functions.storage import inside, read_json, save_json
    from functions.control import run_directory

    from functions.configuration import load_config
    config = load_config(root / args.config)
    if args.cohort:
        config["active_cohort"] = args.cohort
    if args.preflight:
        args.stages = "preflight"
    if args.workers is not None:
        config["workers"] = args.workers
        config["runtime"]["workers"] = args.workers
    if not 1 <= config["workers"] <= 32 or not 1 <= config["threads"] <= 32:
        raise ValueError("工作进程/线程须在 1—32 之间")
    for name in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMBA_NUM_THREADS",
    ]:
        os.environ[name] = str(config["threads"])
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["MPLCONFIGDIR"] = str(root / ".cache/matplotlib")
    os.environ["NUMBA_CACHE_DIR"] = str(root / ".cache/numba")
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    # 二、将唯一标识写进实际命令行，保证停止脚本可复核身份。
    if args.token is None:
        command = [
            sys.executable,
            "-u",
            "-B",
            str(root / "run.py"),
            *sys.argv[1:],
            "--token",
            uuid.uuid4().hex,
        ]
        if args.detach:
            command.remove("--detach")
            directory = run_directory(root, args.run, config["paths"]["results"])
            directory.mkdir(parents=True, exist_ok=True)
            with (directory / "driver.log").open("a") as output:
                process = subprocess.Popen(
                    command,
                    cwd=root,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    start_new_session=True,
                )
            print(
                f'已启动 PID {process.pid}；日志：{directory/"driver.log"}', flush=True
            )
            for iteration in range(50):
                if process.poll() is not None:
                    return process.returncode
                control = directory / "control.json"
                if control.exists() and read_json(control).get("token") == command[-1]:
                    return 0
                time.sleep(0.1)
            return 0
        os.execv(sys.executable, command)
    # 三、实际调度和异常状态记录。
    from functions.runner import run_analysis
    from functions.pipeline import run_pipeline

    try:
        return run_analysis(root, args, config) if args.internal_image else run_pipeline(root, args, config)
    except Exception as exc:
        directory = Path(config.get("output_directory", run_directory(root, args.run, config["paths"]["results"])))
        control = directory / "control.json"
        if control.exists() and read_json(control).get("token") == args.token:
            save_json(
                directory / "progress.json",
                {
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                    "heartbeat": time.time(),
                },
            )
        print(f"运行失败：{type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
