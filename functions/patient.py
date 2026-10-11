"""逐例计算并原子提交；空病灶记录后排除，不生成患者病灶质控图。"""

from __future__ import annotations
import os
import shutil
import time
from pathlib import Path
import numpy as np
from functions.storage import digest, save_json, read_json, verify_files
from functions.reference import load_reference
from functions.space import map_patient, roi_damage
from functions.index import interrupted_counts
from functions.spectral import fixed_support_loss, tissue_loss


def initialize_worker(root: str, threads: int, parent_pid: int) -> None:
    """进程内缓存参考，不在患者之间重读大型索引。"""
    import ctypes
    import signal

    # Linux 父进程异常退出时结束工作进程，防止孤儿进程与续跑同时提交。
    if ctypes.CDLL(None).prctl(1, signal.SIGTERM) != 0 or os.getppid() != parent_pid:
        raise RuntimeError("工作进程无法绑定协调进程的存活状态")
    from threadpoolctl import threadpool_limits
    from numba import set_num_threads

    threadpool_limits(limits=threads)
    set_num_threads(threads)
    process_subject.reference = load_reference(Path(root))


def committed_valid(directory: Path, source_hash: str, signature: str) -> bool:
    if not (directory / "record.json").is_file():
        return False
    try:
        record = read_json(directory / "record.json")
        if record["input_sha256"] != source_hash or record["signature"] != signature:
            return False
        expected = {"result.npz"} if record["automatic_status"] == "passed" else set()
        if (
            record["automatic_status"] not in ["passed", "excluded_empty_lesion"]
            or set(record["outputs"]) != expected
        ):
            return False
        verify_files(directory, record["outputs"])
        return True
    except (OSError, ValueError, KeyError):
        return False


def process_subject(
    root: str, run: str, item: dict, config: dict, signature: str
) -> dict:
    """工作进程不更新全局数据库；只有主控写运行状态。"""
    root = Path(root)
    source = root / item["path"]
    directory = Path(config.get("output_directory", root / "results" / run))
    target = directory / "subjects" / item["subject_id"]
    temporary = directory / "work" / item["subject_id"]
    started = time.monotonic()
    if digest(source) != item["sha256"]:
        raise ValueError("启动后原始病灶发生变化")
    if target.exists():
        if committed_valid(target, item["sha256"], signature):
            return read_json(target / "record.json")
        raise ValueError("已有结果校验失败；请使用新运行名，不覆盖已提交文件")
    if temporary.exists():
        temporary.resolve().relative_to((directory / "work").resolve())
        shutil.rmtree(temporary)  # 仅清理本运行该被试未提交的临时产物
    temporary.mkdir(parents=True)
    reference = process_subject.reference
    empty_file = source.stat().st_size == 0  # 同时记录零字节文件与有效 NIfTI 的空掩膜
    if empty_file:
        mask, summary = None, {"volume_before_mm3": None}
    else:
        mask, summary = map_patient(
            source,
            reference["template"],
            config["maximum_volume_change"],
            allow_empty=True,
            voxel_tolerance=config.get("imaging", {}).get("voxel_tolerance_mm", 1e-4),
        )
    if empty_file or summary["volume_before_mm3"] == 0:
        record = {
            "subject_id": item["subject_id"],
            "input_path": item["path"],
            "input_sha256": item["sha256"],
            "signature": signature,
            "reference_id": config["reference_id"],
            "automatic_status": "excluded_empty_lesion",
            "exclusion_type": "empty_file" if empty_file else "zero_mask",
            "reason": (
                "输入文件为零字节；记录并排除"
                if empty_file
                else "原始病灶掩膜没有非零体素；记录并排除，不生成连接或谐波零值结果"
            ),
            "volume_before_mm3": None if empty_file else 0.0,
            "outputs": {},
            "elapsed_seconds": time.monotonic() - started,
            "lesion_qc": "not_performed_user_request",
        }
        if digest(source) != item["sha256"]:
            raise ValueError("空病灶检查期间输入发生变化")
        save_json(temporary / "record.json", record)
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary.replace(target)
        return record
    hits, rho, anatomical = roi_damage(mask, reference["atlas"], reference["nodes"])
    lost, identifiers, categories = interrupted_counts(
        mask, reference["index"], len(reference["edge_i"])
    )
    ratios, total, _ = fixed_support_loss(lost, reference)
    local, _ = tissue_loss(rho, reference)
    from scipy.sparse.csgraph import connected_components

    remaining = reference["counts"].copy()
    first, second = reference["edge_i"], reference["edge_j"]
    remaining[first, second] -= lost
    remaining[second, first] -= lost
    components = int(
        connected_components(remaining, directed=False, return_labels=False)
    )
    summary.update(anatomical)
    summary.update(
        {
            "subject_id": item["subject_id"],
            "input_path": item["path"],
            "input_sha256": item["sha256"],
            "signature": signature,
            "interrupted_tracks": len(identifiers),
            "paired_interrupted_tracks": int(lost.sum()),
            "unmatched_interrupted_tracks": int(categories[1]),
            "within_roi_interrupted_tracks": int(categories[2]),
            "affected_edges": int(np.count_nonzero(lost)),
            "all_mode_support_summary": total,
            "support_loss_min": float(ratios.min()),
            "support_loss_max": float(ratios.max()),
            "automatic_status": "passed",
            "lesion_qc": "not_performed_user_request",
            "warnings": [],
            "reference_id": config["reference_id"],
            "components_full_135": components,
            "zero_degree_nodes": int(np.sum(remaining.sum(axis=1) == 0)),
        }
    )
    if summary["outside_atlas_voxels"]:
        summary["warnings"].append("病灶含图谱外组织；不能把图谱外覆盖解释为无灰质损伤")
    if len(identifiers) and not lost.sum():
        summary["warnings"].append("命中流线但没有有效脑区端点对；需审核图谱覆盖")
    additional = {"jmax": np.maximum(local, ratios)} if config.get("save_jmax", True) else {}
    np.savez_compressed(
        temporary / "result.npz",
        lost_counts=lost,
        support_loss=ratios,
        w=ratios,
        g=local,
        roi_hit_voxels=hits,
        roi_coverage=rho,
        mode_ids=reference["mode_ids"],
        **additional,
    )
    if digest(source) != item["sha256"]:
        raise ValueError("计算期间原始病灶发生变化；结果未提交")
    summary["elapsed_seconds"] = time.monotonic() - started
    summary["outputs"] = {"result.npz": digest(temporary / "result.npz")}
    for name in summary["outputs"]:
        with (temporary / name).open("rb") as stream:
            os.fsync(stream.fileno())
    save_json(temporary / "record.json", summary)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary.replace(target)
    descriptor = os.open(target.parent, os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return summary
