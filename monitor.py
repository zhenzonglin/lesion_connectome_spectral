#!/usr/bin/env python
"""只读进度监控；显示病例状态和进程存活，不依据静默日志断言卡死。"""
import argparse
import json
import math
import time
from pathlib import Path
from functions.control import run_directory, matching_process
from functions.storage import read_json


def read_status(path: Path) -> dict:
    """首次启动尚未生成状态文件时继续等待，不启动或修改分析。"""
    try:
        return read_json(path)
    except FileNotFoundError:
        return {}


def case_counts(values: dict) -> str:
    return ', '.join(f'{key}={values.get(key, 0)}' for key in
                     ['done', 'failed', 'running', 'queued', 'excluded_empty_lesion'])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--config", default="config.yaml")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--watch", dest="watch", action="store_true", help="持续刷新（默认，兼容旧命令）")
    mode.add_argument("--once", dest="watch", action="store_false", help="只显示一次状态后退出")
    parser.set_defaults(watch=True)
    parser.add_argument("--interval", type=float, help="刷新间隔秒数；默认读取配置（30秒）")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    from functions.configuration import load_config
    config = load_config(root / args.config)
    if args.interval is None:
        args.interval = float(config['runtime'].get('monitor_interval_seconds', 30))
    if not math.isfinite(args.interval) or args.interval < 1:
        raise ValueError("刷新间隔必须为至少 1 秒的有限数值")
    directory = run_directory(root, args.run, config["paths"]["results"])
    while True:
        control = read_status(directory / "control.json")
        progress = read_status(directory / "progress.json")
        progress["process_alive"] = matching_process(root, control) is not None
        progress["heartbeat_age_seconds"] = round(
            time.time() - progress.get("heartbeat", time.time()), 1
        )
        previous_progress = (
            progress.get("heartbeat", 0) < control.get("started", 0)
            or (progress.get("token") is not None and progress["token"] != control.get("token"))
        )
        print(f"\n更新时间 {time.strftime('%H:%M:%S')} | 每 {args.interval:g} 秒刷新", flush=True)
        if not control:
            print("等待总控身份记录；监控不会自行启动分析。", flush=True)
        elif previous_progress or "status" not in progress:
            print(f"等待本次运行刷新进度 | PID {control.get('pid')} | 进程存活 {progress['process_alive']}", flush=True)
        elif 'stages' in progress:
            print(f"阶段进度 | {progress['status']} | 进程存活 {progress['process_alive']}",flush=True)
            for name,stage in progress['stages'].items():
                counts = ''
                stage_path = Path(stage['path']) if stage.get('path') else None
                values = read_status(stage_path / 'progress.json') if stage_path else {}
                cohorts = {}
                if name == 'lesion' and stage_path:
                    for cohort in ['discovery', 'validation']:
                        cohort_values = read_status(stage_path / cohort / 'progress.json')
                        if cohort_values:
                            cohorts[cohort] = cohort_values
                if name == 'lesion' and values:
                    counts=' | '+case_counts(values)
                    counts+=f" | {values.get('subjects_per_minute')} 例/分；预计剩余 {values.get('eta_seconds')} 秒"
                elif values:
                    counts=' | '+', '.join(f'{key}={value}' for key,value in values.items() if key in ['outcome','mode','completed','planned','target','fold','iteration'])
                print(f"{name:14} {stage['status']:12}{counts} {stage.get('reason','')}",flush=True)
                for cohort, values in cohorts.items():
                    print(f"  {cohort:12} {values.get('status', 'unknown'):12} | {case_counts(values)}", flush=True)
        else:
            print(json.dumps(progress, ensure_ascii=False, indent=2), flush=True)
        if not args.watch:
            break
        if control and not progress["process_alive"]:
            print("当前无活动协调进程；继续等待状态更新，按 Ctrl+C 退出监控。", flush=True)
        time.sleep(args.interval)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n监控已退出。")
