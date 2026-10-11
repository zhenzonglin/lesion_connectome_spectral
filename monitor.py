#!/usr/bin/env python
"""只读进度监控；显示病例状态和进程存活，不依据静默日志断言卡死。"""
import argparse
import json
import time
from pathlib import Path
from functions.control import run_directory, matching_process
from functions.storage import read_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    from functions.configuration import load_config
    config = load_config(root / args.config)
    args.interval = args.interval or config['runtime']['monitor_interval_seconds']
    if args.interval < 1:
        raise ValueError("刷新间隔至少 1 秒")
    directory = run_directory(root, args.run, config["paths"]["results"])
    while True:
        progress = read_json(directory / "progress.json")
        control = read_json(directory / "control.json")
        progress["process_alive"] = matching_process(root, control) is not None
        progress["heartbeat_age_seconds"] = round(
            time.time() - progress.get("heartbeat", time.time()), 1
        )
        if 'stages' in progress:
            print(f"阶段进度 | {progress['status']} | 进程存活 {progress['process_alive']}",flush=True)
            for name,stage in progress['stages'].items():
                counts = ''
                status_path = Path(stage.get('path',''))/'progress.json'
                if name == 'lesion' and status_path.is_file():
                    values=read_json(status_path)
                    counts=' | '+', '.join(f'{key}={values.get(key,0)}' for key in ['done','failed','running','queued','excluded_empty_lesion'])
                    counts+=f" | {values.get('subjects_per_minute')} 例/分；预计剩余 {values.get('eta_seconds')} 秒"
                elif status_path.is_file():
                    values=read_json(status_path)
                    counts=' | '+', '.join(f'{key}={value}' for key,value in values.items() if key in ['outcome','mode','completed','planned','target','fold','iteration'])
                print(f"{name:14} {stage['status']:12}{counts} {stage.get('reason','')}",flush=True)
        else:
            print(json.dumps(progress, ensure_ascii=False, indent=2), flush=True)
        if not args.watch or progress["status"] in [
            "completed",
            "completed_no_analyzable_subjects",
            "completed_with_errors",
            "stopped",
            "failed",
            "incomplete",
        ]:
            break
        if not progress["process_alive"]:
            print("协调进程已退出；可检查日志后使用 --resume。")
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
