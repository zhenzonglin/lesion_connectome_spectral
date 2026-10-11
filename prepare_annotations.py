#!/usr/bin/env python
"""注释资源准备/核验。正式运行只读紧凑表；原图和临时构建环境留在 archive。"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to
from scipy.stats import zscore
from functions.storage import digest, fingerprint, read_json, save_json


# 一、分区核查。标签名称、物理坐标和实际体素都检查，不能仅按100行拼表。
def alignment_audit(root: Path) -> tuple[pd.DataFrame, dict]:
    directory = root / "reference/annotation"
    nodes = pd.read_csv(root / "reference/nodes.tsv", sep="\t")
    cortex = nodes.loc[nodes.group == "cortex"].copy()
    labels = pd.read_csv(directory / "raw/labels.txt", sep="\t", header=None)
    source = pd.read_csv(directory / "raw/neurosynth.csv")
    if len(cortex) != 100 or cortex.label.to_list() != labels[1].to_list():
        raise ValueError("项目皮层分区与CBIG正式标签不一致")
    if cortex.label.to_list() != source.iloc[:, 0].to_list():
        raise ValueError("Neurosynth发布表的具名脑区顺序不一致")
    project = nib.load(root / "reference/atlas.nii.gz")
    original = nib.load(directory / "raw/schaefer.nii.gz")
    mapped = resample_from_to(original, project, order=0)
    first, second = np.asarray(project.dataobj), np.asarray(mapped.dataobj)
    rows = []
    for number in cortex["index"]:
        selected, benchmark = first == number, second == number
        shared = int(np.count_nonzero(selected & benchmark))
        rows.append({"index": number, "project_voxels": int(selected.sum()),
                     "original_voxels": int(benchmark.sum()), "intersection_voxels": shared,
                     "dice": 2 * shared / (selected.sum() + benchmark.sum())})
    mismatch = (second > 0) & (first != second)
    unexpected = mismatch & (first <= 100)
    if unexpected.any():
        raise ValueError("存在不能由LQT深部区优先覆盖解释的皮层标签差异")
    audit = pd.DataFrame(rows)
    audit.to_csv(directory / "alignment.tsv", sep="\t", index=False)
    cortex["hemisphere"] = np.where(cortex.label.str.contains("_LH_"), "L", "R")
    cortex[["index", "label", "hemisphere", "x_mm", "y_mm", "z_mm"]].to_csv(
        directory / "cortical_nodes.tsv", sep="\t", index=False)
    summary = {"status": "passed", "matching_label_names": 100,
               "original_compact_tables_accepted": False,
               "overridden_by_deep_parcels_voxels": int(mismatch.sum()),
               "minimum_original_dice": float(audit.dice.min()),
               "mean_original_dice": float(audit.dice.mean()),
               "resolution": "重新从原始影像提取到项目实际LQT皮层体素；不使用原100区表替代"}
    return cortex, summary


def parcel_means(image: nib.spatialimages.SpatialImage,
                 atlas: nib.spatialimages.SpatialImage,
                 remove_background: bool) -> tuple[np.ndarray, np.ndarray]:
    """在原图分辨率最近邻重采样标签；不插值信号。PET删除零背景，NS保留零z值。"""
    image = nib.squeeze_image(image)  # 部分作者PET为 x×y×z×1；只去除末尾单元素维。
    labels = np.asarray(resample_from_to(atlas, image, order=0).dataobj).astype(int)
    values = np.asarray(image.dataobj, dtype=float)
    if values.ndim != 3 or labels.shape != values.shape:
        raise ValueError("注释原图必须为三维且可配准到分区")
    valid = np.isfinite(values) & (labels > 0) & (labels <= 100)
    if remove_background:
        valid &= values != 0
    counts = np.bincount(labels[valid], minlength=101)[1:101]
    sums = np.bincount(labels[valid], weights=values[valid], minlength=101)[1:101]
    if np.any(counts == 0):
        raise ValueError("原始图对部分皮层脑区没有有效体素")
    return sums / counts, counts


# 二、Hansen原图重提取。原作者的聚合运算按实际源码复现。
def build_hansen(root: Path) -> dict:
    archive = root / "archive/annotation_source"
    directory = root / "reference/annotation"
    source = (directory / "raw/make_receptor_matrix.py").read_text()
    names = re.findall(r"scale\+'/([^']+)\.csv'", source)
    urls = read_json(archive / "sources.json")
    atlas = nib.load(root / "reference/atlas.nii.gz")
    extracted, details = [], []
    for name in names:
        files = [path for path in archive.glob(name + ".nii*") if path.is_file()]
        if len(files) != 1:
            extracted.append(np.full(100, np.nan))
            reason = ("OSF公开注册表指向的原图当前返回403/登录页；不使用HTML或原分区均值替代"
                      if name == "VAChT_feobv_hc3_spreng" else "原始图缺失")
            details.append({"source": name, "status": "unavailable", "reason": reason})
            continue
        try:
            image = nib.load(files[0])
            values, coverage = parcel_means(image, atlas, True)
            extracted.append(values)
            details.append({"source": name, "status": "available", "sha256": digest(files[0]),
                            "url": urls.get(files[0].name, ""),
                            "grid": list(image.shape), "affine": image.affine.tolist(),
                            "voxel_size_mm": [float(value) for value in image.header.get_zooms()[:3]],
                            "minimum_parcel_voxels": int(coverage.min())})
        except (nib.filebasedimages.ImageFileError, OSError, ValueError) as exc:
            extracted.append(np.full(100, np.nan))
            details.append({"source": name, "status": "unavailable", "reason": str(exc)})
    data = np.array(extracted).T
    if data.shape != (100, 26):
        raise ValueError("官方Hansen组合清单应含26个来源图")
    source_table = pd.DataFrame(data, columns=names)
    source_table.insert(0, "index", range(1, 101))
    source_table.to_csv(directory / "pet_source_values.tsv", sep="\t", index=False)
    receptors = np.zeros((100, 19))
    receptors[:, 0] = data[:, 0]
    receptors[:, 2:9] = data[:, 3:10]
    receptors[:, 10:14] = data[:, 12:16]
    receptors[:, 15:18] = data[:, 19:22]
    receptors[:, 1] = (zscore(data[:, 1]) * 22 + zscore(data[:, 2]) * 65) / 87
    receptors[:, 9] = (zscore(data[:, 10]) * 37 + zscore(data[:, 11]) * 55) / 92
    receptors[:, 14] = (zscore(data[:, 16]) * 22 + zscore(data[:, 17]) * 28 + zscore(data[:, 18]) * 73) / 123
    # 忠实保留官方实际源码的3、4、1、1；与“按样本数加权”的文字不一致，文档明确标出。
    receptors[:, 18] = (zscore(data[:, 22]) * 3 + zscore(data[:, 23]) * 4 + zscore(data[:, 24]) + zscore(data[:, 25])) / 30
    labels = np.load(directory / "raw/receptor_names.npy", allow_pickle=False).tolist()
    table = pd.DataFrame(receptors, columns=labels)
    table.insert(0, "index", range(1, 101))
    table.to_csv(directory / "hansen.tsv", sep="\t", index=False)
    sensitivity = (zscore(data[:, 22]) * 3 + zscore(data[:, 23]) * 4 + zscore(data[:, 24]) * 5 + zscore(data[:, 25]) * 18) / 30
    pd.DataFrame({"index": range(1, 101), "vacht_sample_weighted": sensitivity}).to_csv(
        directory / "hansen_sensitivity.tsv", sep="\t", index=False)
    available = [label for label, values in zip(labels, receptors.T) if np.isfinite(values).all()]
    record = {"maps_planned": labels, "maps_available": available,
              "sources": details, "background": "finite nonzero PET voxels; native image grid",
              "vacht_formula": "(3*z3+4*z4+z5+z18)/30, as published code",
              "vacht_sensitivity": "(3*z3+4*z4+5*z5+18*z18)/30; not in primary collection"}
    save_json(directory / "hansen_sources.json", record)
    return record


# 三、NS重建只需要一次。兼容别名只改变旧API接口，不改元分析算法。
def neurosynth_runtime(root: Path) -> None:
    runtime = root / "archive/annotation_source/runtime"
    sys.path.insert(0, str(runtime))
    if not (runtime / "neurosynth").is_dir():
        raise RuntimeError("原图重建需要一次性 neurosynth 0.3.8 环境；正式推断无需安装")
    np.round_ = np.round
    np.in1d = np.isin
    nib.spatialimages.SpatialImage.get_data = lambda self: np.asanyarray(self.dataobj)
    nib.spatialimages.SpatialImage.get_header = lambda self: self.header
    nib.spatialimages.SpatialImage.get_affine = lambda self: self.affine
    pd.DataFrame.as_matrix = lambda self: self.to_numpy()
    pd.DataFrame.ix = property(lambda self: self.loc)
    pd.DataFrame.to_sparse = lambda self, *args, **kwargs: self
    pd.DataFrame.to_dense = lambda self: self


def build_neurosynth(root: Path) -> dict:
    """复用原版NS Dataset/MetaAnalysis，固定2018数据库和发布表中的124个术语。"""
    import pickle
    neurosynth_runtime(root)
    from neurosynth import Dataset, MetaAnalysis
    archive = root / "archive/annotation_source"
    directory = root / "reference/annotation"
    terms = pd.read_csv(directory / "raw/neurosynth.csv").columns[1:].to_list()
    cache = archive / "dataset_dense.pkl"
    if cache.exists():
        with cache.open("rb") as stream:
            dataset = pickle.load(stream)
    else:
        image_cache = archive / "dataset_images.pkl"
        if image_cache.exists():
            with image_cache.open("rb") as stream:
                dataset = pickle.load(stream)
        else:
            dataset = Dataset(str(archive / "database/database.txt"))
            with image_cache.open("wb") as stream:
                pickle.dump(dataset, stream, protocol=5)
            print("neurosynth image table cached", flush=True)
        fields = pd.read_csv(archive / "database/features.txt", sep="\t", index_col=0)
        dataset.add_features(fields)
        with cache.open("wb") as stream:
            pickle.dump(dataset, stream, protocol=5)
    atlas = nib.load(root / "reference/atlas.nii.gz")
    original = nib.load(directory / "raw/schaefer.nii.gz")
    extracted, records = [], []
    original_values = []
    raw_maps = archive / "neurosynth_maps"
    raw_maps.mkdir(exist_ok=True)
    for term in terms:
        destination = raw_maps / f"{term}.nii.gz"
        field = term.replace("_", " ")
        selected = dataset.get_studies(features=field)
        if not selected:
            raise ValueError(f"固定术语在原数据库中不存在：{term}")
        if not destination.exists():
            analysis = MetaAnalysis(dataset, selected)
            image = dataset.masker.unmask(analysis.images["association-test_z"], output="image")
            nib.save(image, destination)
        image = nib.load(destination)
        values, coverage = parcel_means(image, atlas, False)
        reference_values, _ = parcel_means(image, original, False)
        original_values.append(reference_values)
        extracted.append(values)
        records.append({"term": term, "selected_studies": len(selected), "sha256": digest(destination),
                        "minimum_parcel_voxels": int(coverage.min())})
        print(f"neurosynth {len(extracted)}/{len(terms)} {term}", flush=True)
    table = pd.DataFrame(np.array(extracted).T, columns=terms)
    table.insert(0, "index", range(1, 101))
    table.to_csv(directory / "neurosynth.tsv", sep="\t", index=False)
    published = pd.read_csv(directory / "raw/neurosynth.csv").iloc[:, 1:].to_numpy()
    correlations = [float(np.corrcoef(first, second)[0, 1])
                    for first, second in zip(np.array(original_values), published.T)]
    record = {"package": "neurosynth 0.3.8", "database_commit": "e8f27c4",
              "database_url": "https://raw.githubusercontent.com/neurosynth/neurosynth-data/e8f27c4/current_data.tar.gz",
              "database_sha256": digest(archive / "current_data.tar.gz"),
              "sphere_radius_mm": 6, "term_threshold": 0.001,
              "map": "association-test_z, unthresholded", "terms": terms,
              "actual_terms": len(terms), "publisher_readme_terms": 123,
              "terms_note": "发布CSV实际含124术语，不任意删列凑成README的123",
              "original_parcellation_correlations_with_published": correlations,
              "sources": records, "background": "zero z values retained"}
    save_json(directory / "neurosynth_sources.json", record)
    return record


def finalize_resources(root: Path) -> dict:
    directory = root / "reference/annotation"
    _, alignment = alignment_audit(root)
    receptors = read_json(directory / "hansen_sources.json")
    cognitive = read_json(directory / "neurosynth_sources.json")
    selected = ["alignment.tsv", "cortical_nodes.tsv", "hansen.tsv", "hansen_sensitivity.tsv", "pet_source_values.tsv",
                "neurosynth.tsv", "hansen_sources.json", "neurosynth_sources.json"]
    selected += [str(path.relative_to(directory)) for path in (directory / "raw").glob("*") if path.is_file()]
    hashes = {name: digest(directory / name) for name in selected}
    manifest = {"resource_id": fingerprint(hashes), "files": hashes, "alignment": alignment,
                "project_nodes_sha256": digest(root / "reference/nodes.tsv"),
                "project_atlas_sha256": digest(root / "reference/atlas.nii.gz"),
                "hansen_available": receptors["maps_available"], "hansen_planned": receptors["maps_planned"],
                "neurosynth_terms": cognitive["actual_terms"],
                "license": "Hansen CC BY-NC-SA 4.0; retain primary source citations",
                "status": "ready" if len(receptors["maps_available"]) == 19 else "partial_resource_availability"}
    save_json(directory / "manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-hansen", action="store_true")
    parser.add_argument("--build-neurosynth", action="store_true")
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.build_hansen:
        result = build_hansen(root)
        print(json.dumps({"hansen_available": len(result["maps_available"]), "planned": 19}))
    if args.build_neurosynth:
        build_neurosynth(root)
    if args.finalize:
        result = finalize_resources(root)
        print(json.dumps({"resource_id": result["resource_id"], "status": result["status"]}))
    if args.verify or not any(vars(args).values()):
        from functions.annotation import load_annotations
        result = load_annotations(root)
        print(json.dumps({"resource_id": result["manifest"]["resource_id"],
                          "hansen_available": len(result["manifest"]["hansen_available"]),
                          "neurosynth_terms": result["manifest"]["neurosynth_terms"]}))


if __name__ == "__main__":
    main()
