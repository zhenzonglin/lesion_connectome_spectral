#!/usr/bin/env python
"""只读进度监控；显示病例状态和进程存活，不依据静默日志断言卡死。"""
import argparse
import json
import math
import time
import unicodedata
from pathlib import Path
from functions.control import run_directory, matching_process
from functions.storage import read_json


def read_status(path: Path) -> dict:
    """首次启动尚未生成状态文件时继续等待，不启动或修改分析。"""
    try:
        return read_json(path)
    except FileNotFoundError:
        return {}


COUNT_KEYS = ('done', 'failed', 'running', 'queued', 'excluded_empty_lesion', 'total')
STATUS_LABELS = {
    'completed': '完成', 'completed_with_errors': '完成但有失败',
    'completed_no_analyzable_subjects': '均被排除', 'running': '运行中',
    'failed': '失败', 'blocked': '受阻', 'pending': '待开始',
    'queued': '待开始', 'preflight': '预检中', 'stopped': '已停止',
    'stopping': '停止中', 'stale': '需更新', 'unknown': '待计数',
    'awaiting_revalidation': '待复核',
}


def patient_counts(values: dict, planned_total: int | None = None) -> dict:
    """未知人数保留为空；排除、失败、运行和排队与成功人数分别记账。"""
    result = dict.fromkeys(COUNT_KEYS)
    result['total'] = planned_total
    ready = all(key in values for key in COUNT_KEYS[:4])
    if ready and (values.get('status') != 'preflight' or 'total' in values):
        result.update({key: int(values.get(key, 0)) for key in COUNT_KEYS[:-1]})
        result['total'] = int(values['total']) if 'total' in values else sum(
            result[key] for key in COUNT_KEYS[:-1])
    return result


def lesion_rows(stage: dict, configured_cohorts: list[str],
                run_started: float | None = None) -> list[tuple]:
    """读取当前病灶阶段两个队列及合计，不扫描或改写原始病灶。"""
    # 续跑开始时仍保留旧路径和计划人数；在阶段重新核查输入前不能展示为当前计数。
    if stage['status'] in ['pending', 'queued', 'stale', 'blocked']:
        return [(name, stage['status'], patient_counts({}))
                for name in ['lesion / 合计', *(f'  {name}' for name in configured_cohorts)]]

    def current_progress(file):
        values = read_status(file)
        if (values and run_started is not None and stage['status'] != 'completed'
                and values.get('heartbeat', 0) < run_started):
            return {'status': 'awaiting_revalidation'}
        return values

    path = Path(stage['path']) if stage.get('path') else None
    totals = stage.get('cohort_totals', {})
    names = list(totals) if totals else configured_cohorts
    cohorts = []
    for name in names:
        folder = path / name if path else None
        values = current_progress(folder / 'progress.json') if folder else {}
        total = totals.get(name)
        fresh = folder is not None and not any(
            (folder / entry).exists() for entry in ['manifest.json', 'state.sqlite', 'subjects'])
        if total is not None and fresh and (not values or values.get('status') == 'preflight'):
            values = dict(values, status=values.get('status', 'pending'),
                          total=total, done=0, failed=0, running=0, queued=total,
                          excluded_empty_lesion=0)
        cohorts.append((f'  {name}', values.get('status', 'unknown'), patient_counts(values, total)))

    # 兼容早期尚未按队列分目录的病灶结果。
    legacy = patient_counts(current_progress(path / 'progress.json')) if path else patient_counts({})
    if not totals and legacy['total'] is not None and all(row[2]['total'] is None for row in cohorts):
        return [('lesion / 合计', stage['status'], legacy)]
    aggregate = {
        key: sum(row[2][key] for row in cohorts)
        if cohorts and all(row[2][key] is not None for row in cohorts) else None
        for key in COUNT_KEYS
    }
    return [('lesion / 合计', stage['status'], aggregate), *cohorts]


def print_stage_table(stages: dict, configured_cohorts: list[str],
                      run_started: float | None = None) -> None:
    rows, details = [], []
    for name, stage in stages.items():
        if name == 'lesion':
            rows.extend(lesion_rows(stage, configured_cohorts, run_started))
        else:
            rows.append((name, stage['status'], patient_counts({})))
            if stage.get('path'):
                values = read_status(Path(stage['path']) / 'progress.json')
                current = ', '.join(f'{key}={value}' for key, value in values.items()
                                    if key in ['outcome', 'mode', 'completed', 'planned',
                                               'target', 'fold', 'iteration'])
                if current:
                    details.append(f'{name}: {current}')
        if stage.get('reason'):
            details.append(f"{name}: {stage['reason']}")

    table = [('阶段/队列', '状态', '已完成', '失败', '进行中', '待完成', '排除', '总人数')]
    table += [(name, STATUS_LABELS.get(status, status),
               *(str(counts[key]) if counts[key] is not None else '—' for key in COUNT_KEYS))
              for name, status, counts in rows]
    def width(value):
        return sum(2 if unicodedata.east_asian_width(char) in 'WF' else 1 for char in value)
    widths = [max(width(row[column]) for row in table) for column in range(8)]
    for index, row in enumerate(table):
        cells = []
        for column, value in enumerate(row):
            padding = ' ' * (widths[column] - width(value))
            cells.append(value + padding if column < 2 else padding + value)
        print(' | '.join(cells), flush=True)
        if index == 0:
            print('-+-'.join('-' * value for value in widths), flush=True)
    print('人数对应病灶逐例任务；待完成为排队人数，失败与排除单列；—表示无逐例人数或尚无计数。', flush=True)
    print('续跑待开始或待复核的队列不沿用上轮计数；人数将在本轮核查后更新。', flush=True)
    for detail in details:
        print(detail, flush=True)


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
        refresh = f'每 {args.interval:g} 秒刷新' if args.watch else '单次查看'
        print(f"\n更新时间 {time.strftime('%H:%M:%S')} | {refresh}", flush=True)
        if not control:
            print("等待总控身份记录；监控不会自行启动分析。", flush=True)
        elif previous_progress or "status" not in progress:
            print(f"等待本次运行刷新进度 | PID {control.get('pid')} | 进程存活 {progress['process_alive']}", flush=True)
        elif 'stages' in progress:
            print(f"阶段进度 | {progress['status']} | 进程存活 {progress['process_alive']} | PID {control.get('pid')}",flush=True)
            cohorts = list(config.get('subjects', {}).get('cohorts', {})) or ['discovery', 'validation']
            print_stage_table(progress['stages'], cohorts, control.get('started'))
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
