"""一、病灶空间核查；二、独立脑区覆盖率；不修改原始病灶。"""

from __future__ import annotations
from pathlib import Path
import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy.ndimage import label


def map_patient(
    path: Path, template, maximum_change: float = 0.01, allow_empty: bool = False,
    voxel_tolerance: float = 1e-4
) -> tuple:
    """按物理坐标最近邻映射；输入须已配准 MNI。本函数不执行解剖配准。"""
    image = nib.load(path)
    data = np.asarray(image.dataobj)
    if data.ndim != 3 or not np.isfinite(data).all() or not np.isin(data, [0, 1]).all():
        raise ValueError("病灶必须是有限数值的三维二值掩膜")
    if (
        not np.isfinite(image.affine).all()
        or abs(np.linalg.det(image.affine[:3, :3])) < 1e-12
    ):
        raise ValueError("病灶仿射矩阵无效")
    if not np.allclose(nib.affines.voxel_sizes(image.affine), 1.0, atol=voxel_tolerance, rtol=0):
        raise ValueError("病灶必须为 1 mm 各向同性体素；请在输入前完成处理")
    same = image.shape == template.shape and np.allclose(
        image.affine, template.affine, atol=1e-5
    )
    mapped = (
        data if same else np.asarray(resample_from_to(image, template, order=0).dataobj)
    )
    mask = mapped > 0
    volume_before = float(data.sum() * abs(np.linalg.det(image.affine[:3, :3])))
    volume_after = float(mask.sum() * abs(np.linalg.det(template.affine[:3, :3])))
    if volume_before == 0 and not allow_empty:
        raise ValueError("真实患者病灶为空；空病灶仅允许自检")
    change = abs(volume_after - volume_before) / volume_before if volume_before else 0.0
    if change > maximum_change:
        raise ValueError(f"重采样体积变化 {change:.2%} 超过 {maximum_change:.2%}")
    coordinates = np.argwhere(mask)
    xyz = nib.affines.apply_affine(template.affine, coordinates)
    components = int(label(mask)[1])
    summary = {
        "volume_before_mm3": volume_before,
        "volume_after_mm3": volume_after,
        "volume_change_fraction": change,
        "resampled": not same,
        "original_shape": list(image.shape),
        "original_affine": image.affine.tolist(),
        "mapped_affine": template.affine.tolist(),
        "components": components,
        "left_voxels": int(np.sum(xyz[:, 0] < 0)),
        "right_voxels": int(np.sum(xyz[:, 0] > 0)),
        "midline_voxels": int(np.sum(xyz[:, 0] == 0)),
        "anatomical_registration": "not_assessed_user_request",
    }
    return mask, summary


def roi_damage(mask: np.ndarray, atlas: np.ndarray, nodes) -> tuple:
    """135 脑区的独立覆盖率；脑干为混合组织，不乘入 SC。"""
    number = len(nodes)
    sizes = np.bincount(atlas.ravel(), minlength=number + 1)[1:]
    hits = np.bincount(atlas[mask], minlength=number + 1)
    rho = np.divide(hits[1:], sizes, out=np.zeros(number), where=sizes > 0)
    brainstem = nodes["group"].eq("brainstem").to_numpy()
    summary = {
        "outside_atlas_voxels": int(hits[0]),
        "covered_atlas_voxels": int(hits[1:].sum()),
        "brainstem_mixed_tissue_voxels": int(hits[1:][brainstem].sum()),
        "gray_parcel_voxels": int(hits[1:][~brainstem].sum()),
        "rho_attenuates_sc": False,
    }
    return hits[1:].astype(np.uint32), rho, summary
