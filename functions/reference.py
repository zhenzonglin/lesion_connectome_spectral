"""冻结参考的加载与数值验收；正式运行不读取 TRK/FIB 或外部工具。"""

from __future__ import annotations
from pathlib import Path
import nibabel as nib
import numpy as np
import pandas as pd
from functions.storage import read_json, verify_files
from functions.index import load_index


def verify_reference(root: Path) -> dict:
    """启动时验证资源全部 SHA256；只读，不重算或修改健康谐波。"""
    record = read_json(root / "reference/manifest.json")
    verify_files(root / "reference", record["files"])
    return record


def load_reference(root: Path, with_index: bool = True) -> dict:
    directory = root / "reference"
    with np.load(directory / "sc.npz") as graph, np.load(
        directory / "modes.npz"
    ) as basis:
        first, second = np.where(np.triu(graph["counts"] > 0, 1))
        reference = {
            name: graph[name].copy() for name in ["weights", "counts", "volumes"]
        }
        reference.update(
            {
                "edge_i": first,
                "edge_j": second,
                "active": basis["active"].copy(),
                "values": basis["values"][basis["positive"]],
                "vectors": basis["vectors"][:, basis["positive"]],
                "mode_ids": np.flatnonzero(basis["positive"]),
                "strength": basis["strength"].copy(),
            }
        )
    reference["nodes"] = pd.read_csv(directory / "nodes.tsv", sep="\t")
    reference["template"] = nib.load(directory / "template.nii.gz")
    atlas_image = nib.load(directory / "atlas.nii.gz")
    if atlas_image.shape != reference["template"].shape or not np.allclose(
        atlas_image.affine, reference["template"].affine
    ):
        raise ValueError("冻结图谱和模板网格不一致")
    reference["atlas"] = np.asarray(atlas_image.dataobj).astype(np.int16)
    if with_index:
        reference["index"] = load_index(directory / "tracks.h5")
        if not np.allclose(reference["index"]["affine"], reference["template"].affine):
            raise ValueError("索引和模板物理坐标不一致")
    return reference


def numerical_checks(reference: dict) -> dict:
    """不删节点、不改权重；零度节点使用已授权的固定排除清单。"""
    graph = reference["weights"]
    if (
        not np.allclose(graph, graph.T, atol=1e-12)
        or np.any(graph < 0)
        or np.any(np.diag(graph))
    ):
        raise ValueError("健康 SC 非对称、负权或存在自环")
    if (
        len(graph),
        int(reference["active"].sum()),
        len(reference["values"]),
        len(reference["edge_i"]),
    ) != (135, 130, 129, 1375):
        raise ValueError("不是本项目冻结的 LQT135 健康参考")
    strength = reference["strength"]
    scale = np.divide(
        1.0, np.sqrt(strength), out=np.zeros_like(strength), where=strength > 0
    )
    laplacian = (
        np.diag(reference["active"].astype(float))
        - scale[:, None] * graph * scale[None, :]
    )
    modes, values = reference["vectors"], reference["values"]
    residual = float(
        np.max(
            np.linalg.norm(laplacian @ modes - modes * values, axis=0)
            / np.maximum(values, 1)
        )
    )
    orthogonal = float(np.max(np.abs(modes.T @ modes - np.eye(len(values)))))
    if max(residual, orthogonal) >= 1e-6:
        raise ValueError("健康谱数值误差超限")
    return {
        "nodes": 135,
        "active_nodes": 130,
        "nonzero_modes": 129,
        "edges": 1375,
        "relative_residual": residual,
        "orthogonality_error": orthogonal,
    }
