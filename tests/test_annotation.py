"""空间注释工程检查；合成病例不会写入正式临床结果。"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import cdist
from functions.annotation import (
    combination_maps, loss_contributions, spatial_weights, moran_surrogates,
    correlation_test, family_bh, load_annotations, run_annotation,
)
from functions.reference import load_reference


@pytest.fixture(scope="module")
def reference():
    return load_reference(Path(__file__).resolve().parents[1], with_index=False)


def test_combination_conservation_and_sign_invariance(reference):
    random = np.random.default_rng(14)
    counts = reference["counts"][reference["edge_i"], reference["edge_j"]]
    lost = np.floor(counts * random.uniform(0, .5, len(counts))).astype(int)
    coverage = random.uniform(0, .2, 135)
    coefficient = random.normal(size=258)
    means, scales = random.uniform(size=(2, 258))
    scales += .2
    result = combination_maps(reference, lost, coverage, coefficient, means, scales)
    assert np.isclose(np.nansum(result["g"]) + np.nansum(result["w"])
                      - result["centering_constant"], result["score"])
    assert np.isnan(result["g"][~reference["active"]]).all()
    flipped = dict(reference, vectors=-reference["vectors"])
    second = combination_maps(flipped, lost, coverage, coefficient, means, scales)
    np.testing.assert_allclose(second["g"], result["g"], equal_nan=True)
    np.testing.assert_allclose(second["w"], result["w"], equal_nan=True)


def test_empty_contributions_and_invalid_input(reference):
    local, network = loss_contributions(reference, np.zeros(1375, dtype=int), np.zeros(135))
    assert np.nansum(local) == 0
    assert np.nansum(network) == 0
    with pytest.raises(ValueError, match="覆盖率"):
        loss_contributions(reference, np.zeros(1375, dtype=int), np.ones(135) * 1.1)


def test_moran_matches_published_singleton_formula():
    # 独立按 BrainSpace moran_randomization 文档中的相关系数公式重建。
    # 本实现用无常数的正交基避免 float32 零特征值误差。
    from scipy.linalg import null_space
    random = np.random.default_rng(88)
    coordinates = random.normal(size=(10, 3))
    hemisphere = np.array(["L"] * 5 + ["R"] * 5)
    values = random.normal(size=10)
    weights = spatial_weights(coordinates, hemisphere)
    nulls, _ = moran_surrogates(values, weights, hemisphere, n_perm=30, seed=101)
    independent_random = np.random.default_rng(101)
    expected = np.empty_like(nulls)
    for label in ("L", "R"):
        selected = np.flatnonzero(hemisphere == label)
        block = weights[np.ix_(selected, selected)]
        basis = null_space(np.ones((1, 5)))
        _, vectors = np.linalg.eigh(basis.T @ block @ basis)
        mem = basis @ vectors
        correlations = 1 - cdist(values[selected][None, :], mem.T, metric="correlation").ravel()
        signs = independent_random.choice(np.array([-1., 1.]), size=(4, 30))
        expected[selected] = values[selected].mean() + mem @ (correlations[:, None] * signs) * (np.sqrt(4) * values[selected].std(ddof=1))
    np.testing.assert_allclose(nulls, expected, atol=1e-12)


def test_spatial_null_preserves_hemisphere_variance_and_moran():
    random = np.random.default_rng(2)
    coordinates = random.normal(size=(12, 3))
    hemisphere = np.repeat(["L", "R"], 6)
    values = random.normal(size=12)
    weights = spatial_weights(coordinates, hemisphere)
    nulls, _ = moran_surrogates(values, weights, hemisphere, 100, 3)
    for label in ("L", "R"):
        selected = hemisphere == label
        block = weights[np.ix_(selected, selected)]
        target = values[selected] - values[selected].mean()
        centered = nulls[selected] - nulls[selected].mean(axis=0)
        np.testing.assert_allclose(nulls[selected].mean(axis=0), values[selected].mean(), atol=1e-12)
        np.testing.assert_allclose(np.sum(centered ** 2, axis=0), target @ target, atol=1e-12)
        np.testing.assert_allclose(np.sum(centered * (block @ centered), axis=0), target @ block @ target, atol=1e-12)


def test_spatial_pvalue_and_failed_family_denominator():
    values = np.arange(10, dtype=float)
    nulls = np.c_[values, -values]
    correlations, probabilities = correlation_test(values, values[:, None], nulls)
    np.testing.assert_allclose(correlations, [1])
    np.testing.assert_allclose(probabilities, [1])
    np.testing.assert_allclose(family_bh(np.array([.01, np.nan, .04])), [.03, np.nan, .06], equal_nan=True)


def test_annotation_resource_alignment_and_gating(tmp_path):
    root = Path(__file__).resolve().parents[1]
    resources = load_annotations(root)
    assert len(resources["nodes"]) == 100
    assert resources["manifest"]["alignment"]["original_compact_tables_accepted"] is False
    pd.DataFrame({"index": range(1, 136), "value": np.sin(np.arange(135))}).to_csv(
        tmp_path / "regions.tsv", sep="\t", index=False)
    (tmp_path / "maps_manifest.json").write_text(json.dumps({"maps": [
        {"pattern_id": "synthetic", "channel": "g", "map_file": "regions.tsv", "validated": True},
        {"pattern_id": "unreplicated", "channel": "w", "map_file": "regions.tsv", "validated": False},
    ]}))
    summary = run_annotation(root, {"annotation": {"permutations": 49}}, tmp_path, tmp_path / "output")
    assert summary["n_perm"] == 49
    table = pd.read_csv(tmp_path / "output/associations.tsv", sep="\t")
    assert table.loc[table.pattern_id == "unreplicated", "p"].isna().all()
    assert set(table.family) == {"neurosynth", "hansen"}
    assert len(table) == 2 * ((resources["hansen"].shape[1] - 1) + (resources["neurosynth"].shape[1] - 1))
