"""一、精确线段遍历；二、体素倒排索引；三、唯一流线的脑区对计数。"""

from __future__ import annotations
from pathlib import Path
import h5py
import numpy as np
from numba import njit


@njit(cache=True)
def append_segment(start, end, shape, buffer, position, track):
    """复用试运行 segment_hit 的 DDA 边界规则；写出体素与流线组合键。"""
    delta = end - start  # 体素中心坐标之差；边界在半整数处
    lower, upper = 0.0, 1.0
    for axis in range(3):
        if abs(delta[axis]) < 1e-14:
            if start[axis] < -0.5 or start[axis] >= shape[axis] - 0.5:
                return position
        else:
            first = (-0.5 - start[axis]) / delta[axis]
            last = (shape[axis] - 0.5 - start[axis]) / delta[axis]
            if first > last:
                first, last = last, first
            lower, upper = max(lower, first), min(upper, last)
    if lower > upper:
        return position
    cell = np.floor(start + delta * (lower + 1e-12 * (upper - lower)) + 0.5).astype(
        np.int64
    )
    step = np.zeros(3, dtype=np.int64)
    crossing, interval = np.full(3, np.inf), np.full(3, np.inf)
    for axis in range(3):
        cell[axis] = min(max(cell[axis], 0), shape[axis] - 1)
        if delta[axis] > 1e-14:
            step[axis] = 1
            crossing[axis] = (cell[axis] + 0.5 - start[axis]) / delta[axis]
            interval[axis] = 1.0 / delta[axis]
        elif delta[axis] < -1e-14:
            step[axis] = -1
            crossing[axis] = (cell[axis] - 0.5 - start[axis]) / delta[axis]
            interval[axis] = -1.0 / delta[axis]
    for iteration in range(sum(shape) + 6):
        voxel = (cell[0] * shape[1] + cell[1]) * shape[2] + cell[2]
        buffer[position] = (np.uint64(voxel) << np.uint64(32)) | np.uint64(track)
        position += 1
        next_value = min(crossing[0], crossing[1], crossing[2])
        if next_value > upper + 1e-12 or not np.isfinite(next_value):
            return position
        for axis in range(3):
            if crossing[axis] <= next_value + 1e-12:
                cell[axis] += step[axis]
                crossing[axis] += interval[axis]
        if np.any(cell < 0) or np.any(cell >= shape):
            return position
    return position


@njit(cache=True)
def collect_memberships(points, offsets, identifiers, transform, shape):
    """批量生成倒排键；点坐标为原生 DSI 体素坐标，短流线不参与。"""
    mapped = np.empty(points.shape, dtype=np.float64)
    for point in range(len(points)):
        for axis in range(3):
            mapped[point, axis] = transform[axis, 3]
            for coordinate in range(3):
                mapped[point, axis] += (
                    transform[axis, coordinate] * points[point, coordinate]
                )
    capacity = 0  # 每条线段的保守穿越体素数上限
    for track in range(len(identifiers)):
        for point in range(offsets[track], offsets[track + 1] - 1):
            capacity += 4 + int(
                np.ceil(np.abs(mapped[point + 1] - mapped[point])).sum()
            )
    buffer = np.empty(capacity, dtype=np.uint64)
    position = 0
    for track in range(len(identifiers)):
        for point in range(offsets[track], offsets[track + 1] - 1):
            position = append_segment(
                mapped[point],
                mapped[point + 1],
                shape,
                buffer,
                position,
                identifiers[track],
            )
    return buffer[:position]


@njit(cache=True)
def query_ids(voxels, pointers, tracks, total):
    """所有病灶体素命中的流线取并集，每条只计一次。"""
    selected = np.zeros(total, dtype=np.bool_)
    for voxel in voxels:
        for position in range(pointers[voxel], pointers[voxel + 1]):
            selected[tracks[position]] = True
    return np.flatnonzero(selected)


def load_index(path: Path) -> dict:
    """每个工作进程只加载一次，避免逐体素 HDF5 随机解压。"""
    with h5py.File(path, "r") as source:
        return {
            name: source[name][:]
            for name in [
                "pointers",
                "tracks",
                "edge_ids",
                "endpoints",
                "classes",
                "shape",
                "affine",
            ]
        }


def interrupted_counts(mask: np.ndarray, index: dict, edges: int) -> tuple:
    if tuple(mask.shape) != tuple(index["shape"]):
        raise ValueError("病灶与倒排索引尺寸不一致")
    identifiers = query_ids(
        np.flatnonzero(mask.ravel()),
        index["pointers"],
        index["tracks"],
        len(index["edge_ids"]),
    )
    paired = index["edge_ids"][identifiers]
    counts = np.bincount(paired[paired >= 0], minlength=edges).astype(np.uint32)
    classes = np.bincount(index["classes"][identifiers], minlength=4)
    return counts, identifiers, classes
