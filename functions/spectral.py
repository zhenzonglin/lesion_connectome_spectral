"""一、计数到加权 SC；二、固定健康模态支撑损失；三、按需节点贡献。"""

from __future__ import annotations
import numpy as np


def patient_weights(healthy_counts, lost_counts, volumes):
    """所有矩阵沿用固定健康节点；脑区体积单位为 mm³。"""
    if np.any(lost_counts < 0) or np.any(lost_counts > healthy_counts):
        raise ValueError("中断计数必须位于零和健康计数之间")
    return (healthy_counts - lost_counts) / (volumes[:, None] + volumes[None, :])


def fixed_support_loss(lost, reference, contributions: bool = False) -> tuple:
    """边形式与原二次型严格等价；需要时才生成脑区×模态贡献。"""
    first, second = reference["edge_i"], reference["edge_j"]
    counts = reference["counts"][first, second]
    if lost.shape != counts.shape or np.any(lost < 0) or np.any(lost > counts):
        raise ValueError("患者中断计数不满足健康连接上界")
    delta = lost.astype(np.float64) / (
        reference["volumes"][first] + reference["volumes"][second]
    )
    values = reference["values"]  # 129 个非零特征值
    modes = reference["vectors"]  # 135×129；仅健康有连接节点参加计算
    strength = reference["strength"]
    scale = np.divide(
        1.0, np.sqrt(strength), out=np.zeros_like(strength), where=strength > 0
    )
    difference = modes[first] * scale[first, None] - modes[second] * scale[second, None]
    energy = delta[:, None] * difference**2
    ratios = energy.sum(axis=0) / values
    if np.min(ratios) < -1e-8 or np.max(ratios) > 1 + 1e-8:
        raise ValueError("谐波支撑损失超出 [0,1]")
    spatial = None
    if contributions:
        spatial = np.zeros_like(modes)
        np.add.at(spatial, first, energy / (2 * values))
        np.add.at(spatial, second, energy / (2 * values))
        if not np.allclose(spatial.sum(axis=0), ratios, atol=1e-10, rtol=1e-8):
            raise ValueError("脑区贡献与模态损失不守恒")
        spatial[~reference["active"]] = np.nan
    return ratios, float(energy.sum() / values.sum()), spatial


def tissue_loss(rho: np.ndarray, reference: dict, contributions: bool = False) -> tuple:
    """局部组织损伤 G；rho 为135脑区覆盖率，平方载荷消除模态符号任意性。"""
    rho = np.asarray(rho, dtype=float)
    modes = reference["vectors"]  # 135×129 的冻结健康模态
    if rho.shape != (modes.shape[0],) or not np.isfinite(rho).all():
        raise ValueError("脑区覆盖率维度或有限性错误")
    if np.any((rho < 0) | (rho > 1)):
        raise ValueError("脑区覆盖率须位于 [0,1]")
    spatial = rho[:, None] * modes**2
    loss = spatial.sum(axis=0)
    if np.any((loss < -1e-10) | (loss > 1+1e-10)):
        raise ValueError("G 超出 [0,1]，请检查健康模态归一化")
    if contributions:
        spatial[~reference["active"]] = np.nan
    return loss, spatial if contributions else None
