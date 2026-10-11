"""病例发现、运行身份和临床特征接口；不读取或拟合临床结局。"""

from __future__ import annotations
import csv
import fcntl
import importlib.metadata
import re
from pathlib import Path
import numpy as np
import pandas as pd
from functions.storage import digest, fingerprint, inside, read_json, save_json


def configured_cohorts(config: dict) -> list[str]:
    """返回配置实际启用的队列；可由 --cohort 限定单个队列。"""
    paths = config.get("paths", {}).get("lesion", "lesion")
    subjects = config.get("subjects", {})
    available = ["discovery", "validation"] if isinstance(paths, dict) or subjects.get("cohorts") else ["discovery"]
    active = config.get("active_cohort")
    if active is not None:
        if active not in available:
            raise ValueError(f"未知 active_cohort：{active}")
        return [active]
    return available


def cohort_settings(config: dict, cohort: str = "discovery") -> tuple[Path | str, dict]:
    """返回队列病灶目录设置；兼容旧版单目录配置。"""
    paths = config.get("paths", {})
    configured_path = paths.get("lesion", "lesion")
    if isinstance(configured_path, dict):
        lesion_path = configured_path.get(cohort)
        if not lesion_path:
            raise ValueError(f"paths.lesion 缺少 {cohort} 路径")
    else:
        cohort_map = config.get("subjects", {}).get("cohorts", {})
        lesion_path = f"{configured_path}/{cohort}" if cohort in cohort_map else configured_path
    subjects = config.get("subjects", {})
    settings = dict(subjects)
    settings.update(subjects.get("cohorts", {}).get(cohort, {}))
    settings.setdefault("pattern", "sub-*_lesion.nii.gz")
    return lesion_path, settings


def extract_subject_id(filename: str, settings: dict) -> str:
    """按正则或前后缀提取 ID；正则必须含 subject_id 命名组。"""
    rule = settings.get("id_rule", {}) or {}
    regex = settings.get("id_regex") or rule.get("id_regex")
    if regex:
        match = re.fullmatch(regex, filename)
        if not match or "subject_id" not in match.groupdict():
            raise ValueError(f"病灶文件名不符合 id_regex：{filename}")
        return match.group("subject_id")
    pattern = settings.get("pattern", "sub-*_lesion.nii.gz")
    prefix = rule.get("prefix", pattern.split("*", 1)[0] if "*" in pattern else "")
    suffix = rule.get("suffix", pattern.split("*", 1)[1] if "*" in pattern else "")
    if pattern.count("*") != 1 or any(symbol in pattern for symbol in "?[]/\\"):
        raise ValueError("病灶文件模式必须包含一个 *，或使用含 subject_id 的 id_regex")
    if not filename.startswith(prefix) or (suffix and not filename.endswith(suffix)):
        raise ValueError(f"病灶文件名不符合 pattern：{filename}")
    stop = len(filename) - len(suffix) if suffix else None
    identifier = filename[len(prefix):stop]
    if not identifier:
        raise ValueError(f"病灶文件名提取到空 subject_id：{filename}")
    return identifier


def discover(root: Path, subjects: str, config: dict, cohort: str = "discovery") -> list[dict]:
    """扫描一个队列的病灶；不改变大小写、前导零或 ses 文本。"""
    items = []
    configured_path, options = cohort_settings(config, cohort)
    lesion_dir = (root / configured_path).resolve()
    if not lesion_dir.is_dir():
        raise FileNotFoundError(f"{cohort} 病灶目录不存在：{lesion_dir}")
    participants = options.get("participants")
    if participants:
        records = pd.read_csv(root / participants, sep="\t", dtype=str, keep_default_na=False)
        if not {"subject_id", "lesion"}.issubset(records.columns):
            raise ValueError("participants.tsv 必须包含 subject_id、lesion")
        candidates = [(row.subject_id, (root / row.lesion).resolve()) for row in records.itertuples()]
    else:
        pattern = options.get("pattern", "sub-*_lesion.nii.gz")
        candidates = []
        if (pattern.count("*") != 1 or any(symbol in pattern for symbol in "?[]/\\")) and not options.get("id_regex") and not options.get("id_rule", {}).get("id_regex"):
            raise ValueError("病灶文件模式必须包含恰好一个 *，或配置 id_regex")
        has_regex = bool(options.get("id_regex") or options.get("id_rule", {}).get("id_regex"))
        # 正则负责完整文件名匹配，不能先用默认的 sub-*_lesion.nii.gz 把复杂文件名漏掉。
        candidate_paths = lesion_dir.glob("*") if has_regex else lesion_dir.glob(pattern if "*" in pattern else "*")
        for path in sorted(candidate_paths):
            if not path.is_file():
                continue
            try:
                identifier = extract_subject_id(path.name, options)
            except ValueError:
                if options.get("id_regex") or options.get("id_rule", {}).get("id_regex"):
                    continue
                raise
            candidates.append((identifier, path))
    selected = (
        None
        if subjects.lower() == "all"
        else set(subjects.split(","))
    )
    if selected is not None and len(selected) != len(subjects.split(",")):
        raise ValueError("选择清单包含重复 ID")
    seen = set()
    for subject, path in candidates:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", subject):
            raise ValueError("患者编号含不安全字符")
        if subject in seen:
            raise ValueError(
                f"同一患者有多个病灶文件，请通过 participants.tsv 区分：{subject}"
            )
        seen.add(subject)
        if selected is not None and subject not in selected:
            continue
        path = path.resolve()
        path.relative_to(lesion_dir)
        if not path.is_file():
            raise FileNotFoundError(path)
        if config.get("space") != "MNI152NLin6Asym":
            raise ValueError("输入空间必须与冻结参考 MNI152NLin6Asym 一致")
        if not config.get("imaging", {}).get("space_confirmed", False) and config["space"] not in path.name:
            raise ValueError("未确认病灶 MNI 空间：请在 imaging.space_confirmed 显式确认配准来源")
        items.append(
            {
                "cohort": cohort,
                "subject_id": subject,
                "path": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
                "sha256": digest(path),
            }
        )
    if selected is not None and selected - seen:
        raise ValueError(f"缺少病灶：{sorted(selected-seen)}")
    if not items:
        raise ValueError(f"{cohort} 队列没有可分析的病灶")
    return items


def method_identity(root: Path, config: dict, reference: dict) -> dict:
    names = ["space.py", "index.py", "spectral.py", "reference.py", "patient.py"]
    versions = {
        name: importlib.metadata.version(name)
        for name in ["numpy", "scipy", "nibabel", "numba", "h5py"]
    }
    return {
        "reference_id": reference["reference_id"],
        "parameters": {key: config.get(key) for key in ["space", "maximum_volume_change", "save_jmax", "imaging"]},
        "core_code": {name: digest(root / "functions" / name) for name in names},
        "versions": versions,
    }


def export_tables(root: Path, directory: Path, jobs: list) -> None:
    """只汇总已成功提交的病例；失败项单列，不用零代替。"""
    spectra, coverage, summaries, excluded = [], [], [], []
    for subject, status, error in jobs:
        if status == "excluded_empty_lesion":
            record = read_json(directory / "subjects" / subject / "record.json")
            excluded.append(
                {
                    key: record[key]
                    for key in [
                        "subject_id",
                        "input_path",
                        "input_sha256",
                        "automatic_status",
                        "reason",
                    ]
                }
            )
        if status != "done":
            continue
        target = directory / "subjects" / subject
        record = read_json(target / "record.json")
        with np.load(target / "result.npz") as data:
            spectra.append(
                {
                    "subject_id": subject,
                    **{
                        f"w_{int(k):03d}": float(v)
                        for k, v in zip(data["mode_ids"], data["support_loss"])
                    },
                    **{f"g_{int(k):03d}": float(v) for k,v in zip(data["mode_ids"],data["g"])},
                    **({f"jmax_{int(k):03d}": float(v) for k,v in zip(data["mode_ids"],data["jmax"])} if "jmax" in data else {}),
                }
            )
            coverage.append(
                {
                    "subject_id": subject,
                    **{
                        f"roi_{i+1}": float(v)
                        for i, v in enumerate(data["roi_coverage"])
                    },
                }
            )
        summaries.append(
            {
                key: record[key]
                for key in [
                    "subject_id",
                    "volume_after_mm3",
                    "interrupted_tracks",
                    "paired_interrupted_tracks",
                    "affected_edges",
                    "all_mode_support_summary",
                    "outside_atlas_voxels",
                    "elapsed_seconds",
                ]
            }
        )
        summaries[-1]["lesion_volume_mm3"] = record["volume_before_mm3"]
    with np.load(root / "reference/modes.npz") as basis:
        modes = np.flatnonzero(basis["positive"])
    for name, rows, columns in [
        ("mode_features", spectra, ["subject_id"] + [f"{channel}_{k:03d}" for channel in ["g","w"] for k in modes] + ([f"jmax_{k:03d}" for k in modes] if spectra and "jmax_001" in spectra[0] else [])),
        (
            "roi_coverage",
            coverage,
            ["subject_id"] + [f"roi_{i}" for i in range(1, 136)],
        ),
    ]:
        pd.DataFrame(rows, columns=columns).to_csv(
            directory / (name + ".tsv.gz"), sep="\t", index=False, compression="gzip"
        )
    summary_columns = [
        "subject_id",
        "lesion_volume_mm3",
        "volume_after_mm3",
        "interrupted_tracks",
        "paired_interrupted_tracks",
        "affected_edges",
        "all_mode_support_summary",
        "outside_atlas_voxels",
        "elapsed_seconds",
    ]
    pd.DataFrame(summaries, columns=summary_columns).to_csv(
        directory / "summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        excluded,
        columns=[
            "subject_id",
            "input_path",
            "input_sha256",
            "automatic_status",
            "reason",
        ],
    ).to_csv(directory / "excluded.tsv", sep="\t", index=False)
    pd.DataFrame(jobs, columns=["subject_id", "calculation_status", "error"]).to_csv(
        directory / "status.tsv", sep="\t", index=False
    )
    # 数值和临床以 subject_id 连接；本阶段没有临床数据进入计算。
