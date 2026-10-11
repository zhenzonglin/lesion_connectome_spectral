"""独立数值与几何测试：不调用外部工具，不需要真实患者数据。"""

from pathlib import Path
import numpy as np
import nibabel as nib
import pytest
from tests.oracle import segment_hit
from functions.index import collect_memberships
from functions.reference import load_reference, numerical_checks
from functions.spectral import fixed_support_loss
from functions.space import map_patient, roi_damage


def cells_for_track(points, shape):
    keys = collect_memberships(
        np.asarray(points, dtype=float),
        np.array([0, len(points)]),
        np.array([0]),
        np.eye(4),
        np.array(shape),
    )
    return np.unique((keys >> np.uint64(32)).astype(np.int64))


def test_sparse_small_and_repeated_path():
    shape = (8, 9, 10)
    points = [[0, 4, 5], [7, 4, 5], [0, 4, 5]]
    cells = cells_for_track(points, shape)
    assert np.ravel_multi_index((3, 4, 5), shape) in cells
    assert len(cells) == 8


def test_random_and_boundary_segments_against_original_oracle():
    random = np.random.default_rng(17)
    shape = (8, 9, 10)
    examples = [
        np.array([[0, 0, 0], [7, 8, 9]]),
        np.array([[-0.5, 2, 2], [7.5, 2, 2]]),
        np.array([[2.5, 2.5, 2.5], [2.5, 2.5, 2.5]]),
        np.array([[-5, -5, -5], [-1, -1, -1]]),
    ]
    examples += [random.uniform(-2, 11, (2, 3)) for _ in range(80)]
    for points in examples:
        for oriented in [points, points[::-1].copy()]:
            cells = cells_for_track(oriented, shape)
            for iteration in range(15):
                mask = random.random(shape) < 0.015
                expected = bool(segment_hit(oriented[0], oriented[1], mask))
                assert bool(mask.ravel()[cells].any()) == expected


def test_fixed_scale_matches_independent_matrix_formula():
    reference = load_reference(Path(__file__).resolve().parents[1], with_index=False)
    numerical_checks(reference)
    first, second = reference["edge_i"], reference["edge_j"]
    healthy = reference["counts"][first, second]
    for fraction in [0.0, 0.25, 1.0]:
        ratios, _, spatial = fixed_support_loss(healthy * fraction, reference, True)
        assert np.allclose(ratios, fraction, atol=1e-10)
        assert np.allclose(np.nansum(spatial, axis=0), ratios, atol=1e-10)
        assert np.isnan(spatial[~reference["active"]]).all()
    random = np.random.default_rng(13)
    lost = np.floor(healthy * random.random(len(healthy)))
    ratios, _, _ = fixed_support_loss(lost, reference)
    delta = np.zeros_like(reference["weights"])
    delta[first, second] = lost / (
        reference["volumes"][first] + reference["volumes"][second]
    )
    delta += delta.T
    strength = reference["strength"]
    scale = np.divide(
        1, np.sqrt(strength), out=np.zeros_like(strength), where=strength > 0
    )
    operator = scale[:, None] * (np.diag(delta.sum(axis=1)) - delta) * scale[None, :]
    expected = (
        np.diag(reference["vectors"].T @ operator @ reference["vectors"])
        / reference["values"]
    )
    assert np.allclose(ratios, expected, atol=1e-10)


def test_space_rejects_wrong_or_empty_input(tmp_path):
    template = nib.Nifti1Image(np.zeros((5, 5, 5)), np.eye(4))
    path = tmp_path / "mask.nii.gz"
    nib.save(template, path)
    with pytest.raises(ValueError, match="为空"):
        map_patient(path, template)
    data = np.zeros((5, 5, 5))
    data[2, 2, 2] = 1
    affine = np.eye(4)
    affine[0, 3] = 100
    nib.save(nib.Nifti1Image(data, affine), path)
    with pytest.raises(ValueError, match="体积变化"):
        map_patient(path, template)
    data[1, 1, 1] = 0.5
    nib.save(nib.Nifti1Image(data, np.eye(4)), path)
    with pytest.raises(ValueError, match="二值"):
        map_patient(path, template)


@pytest.mark.parametrize("kind", ["empty_file", "zero_mask"])
def test_empty_patient_is_recorded_and_excluded(tmp_path, monkeypatch, kind):
    """真实患者入口排除空输入；不会访问流线索引或生成零值谱。"""
    from functions.patient import process_subject, committed_valid
    from functions.storage import digest, read_json

    template = nib.Nifti1Image(np.zeros((5, 5, 5), dtype=np.uint8), np.eye(4))
    source = tmp_path / "lesion" / "test.nii.gz"
    source.parent.mkdir()
    if kind == "empty_file":
        source.write_bytes(b"")
    else:
        nib.save(template, source)
    monkeypatch.setattr(
        process_subject, "reference", {"template": template}, raising=False
    )
    item = {
        "subject_id": "TEST",
        "path": "lesion/test.nii.gz",
        "sha256": digest(source),
    }
    record = process_subject(
        str(tmp_path),
        "empty",
        item,
        {"maximum_volume_change": 0.01, "reference_id": "unit"},
        "unit_signature",
    )
    directory = tmp_path / "results/empty/subjects/TEST"
    assert (
        record["automatic_status"] == "excluded_empty_lesion"
        and record["exclusion_type"] == kind
    )
    assert {path.name for path in directory.iterdir()} == {"record.json"}
    assert committed_valid(directory, item["sha256"], "unit_signature")
