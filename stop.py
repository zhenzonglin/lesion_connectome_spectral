#!/usr/bin/env python
"""安全中断：默认当前病例完成后停止；--now 立即停止自有工作进程。"""
import argparse
from pathlib import Path
from functions.control import request_stop


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--now", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    from functions.configuration import load_config
    config = load_config(root / args.config)
    print(request_stop(root, args.run, args.now, config["paths"]["results"]))


if __name__ == "__main__":
    main()
