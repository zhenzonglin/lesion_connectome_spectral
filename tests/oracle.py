"""试运行已验证的独立完整线段相交实现，仅用于回归验收。"""

import numpy as np
from numba import njit


@njit(cache=True)
def segment_hit(start: np.ndarray, end: np.ndarray, mask: np.ndarray) -> bool:
    """遍历线段穿过的全部体素。

    输入为体素中心坐标（整数对应中心）；立方体边界位于半整数。
    使用三维 DDA，不以固定步长近似整条路径。
    """
    delta = end - start  # 两端差，体素单位
    lower, upper = 0.0, 1.0  # 线段进入/离开影像范围的参数
    for axis in range(3):
        if abs(delta[axis]) < 1e-14:
            if start[axis] < -0.5 or start[axis] >= mask.shape[axis] - 0.5:
                return False
        else:
            first = (-0.5 - start[axis]) / delta[axis]
            last = (mask.shape[axis] - 0.5 - start[axis]) / delta[axis]
            if first > last:
                first, last = last, first
            lower = max(lower, first)
            upper = min(upper, last)
    if lower > upper:
        return False
    position = start + delta * (lower + 1e-12 * (upper - lower))
    cell = np.floor(position + 0.5).astype(np.int64)
    step = np.zeros(3, dtype=np.int64)
    next_crossing = np.full(3, np.inf)
    interval = np.full(3, np.inf)
    for axis in range(3):
        cell[axis] = min(max(cell[axis], 0), mask.shape[axis] - 1)
        if delta[axis] > 1e-14:
            step[axis] = 1
            next_crossing[axis] = (cell[axis] + 0.5 - start[axis]) / delta[axis]
            interval[axis] = 1.0 / delta[axis]
        elif delta[axis] < -1e-14:
            step[axis] = -1
            next_crossing[axis] = (cell[axis] - 0.5 - start[axis]) / delta[axis]
            interval[axis] = -1.0 / delta[axis]
    for iteration in range(sum(mask.shape) + 6):
        if mask[cell[0], cell[1], cell[2]]:
            return True
        next_value = min(next_crossing[0], next_crossing[1], next_crossing[2])
        if next_value > upper + 1e-12 or not np.isfinite(next_value):
            return False
        for axis in range(3):
            if next_crossing[axis] <= next_value + 1e-12:
                cell[axis] += step[axis]
                next_crossing[axis] += interval[axis]
        if (
            np.any(cell < 0)
            or cell[0] >= mask.shape[0]
            or cell[1] >= mask.shape[1]
            or cell[2] >= mask.shape[2]
        ):
            return False
    return False
