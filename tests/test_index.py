"""冻结健康索引的模拟阳性对照和灰白质结果分离检查。"""

from pathlib import Path
import numpy as np
import pytest
from functions.reference import load_reference
from functions.space import map_patient, roi_damage
from functions.index import interrupted_counts
from functions.spectral import fixed_support_loss


@pytest.fixture(scope="module")
def reference():
    return load_reference(Path(__file__).resolve().parents[1])


@pytest.mark.parametrize(
    "name,total,paired", [("pons", 556, 550), ("cerebellum", 497, 379)]
)
def test_archived_infratentorial_controls(reference, name, total, paired):
    directory = Path(__file__).parent / "fixtures" / name
    mask, _ = map_patient(directory / "lesion.nii.gz", reference["template"])
    counts, identifiers, _ = interrupted_counts(
        mask, reference["index"], len(reference["edge_i"])
    )
    assert np.array_equal(identifiers, np.load(directory / "expected_ids.npy"))
    assert len(identifiers) == total and counts.sum() == paired
    ratios, _, _ = fixed_support_loss(counts, reference)
    assert ratios.max() > 0


def test_empty_and_gray_white_separation(reference):
    mask = np.zeros(reference["atlas"].shape, dtype=bool)
    counts, identifiers, _ = interrupted_counts(
        mask, reference["index"], len(reference["edge_i"])
    )
    ratios, _, _ = fixed_support_loss(counts, reference)
    assert not len(identifiers) and np.max(np.abs(ratios)) <= 1e-8
    pointer = reference["index"]["pointers"]
    available = np.flatnonzero(np.diff(pointer) > 0)
    selected = None
    for voxel in available[reference["atlas"].ravel()[available] == 0]:
        ids = reference["index"]["tracks"][pointer[voxel] : pointer[voxel + 1]]
        if np.any(reference["index"]["edge_ids"][ids] >= 0):
            selected = voxel
            break
    assert selected is not None
    mask.ravel()[selected] = True
    counts, _, _ = interrupted_counts(
        mask, reference["index"], len(reference["edge_i"])
    )
    _, rho, _ = roi_damage(mask, reference["atlas"], reference["nodes"])
    assert counts.sum() > 0 and rho.sum() == 0
    mask[:] = False
    mask.ravel()[np.flatnonzero(reference["atlas"].ravel() > 0)[0]] = True
    _, rho, _ = roi_damage(mask, reference["atlas"], reference["nodes"])
    assert rho.sum() > 0  # 覆盖率单独计算，不作为谱函数的输入。
