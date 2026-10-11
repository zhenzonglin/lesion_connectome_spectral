"""空间对象、组合贡献及皮层空间零模型；不参与临床权重选择。"""

from __future__ import annotations

from pathlib import Path
from collections.abc import Iterable
import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.spatial.distance import cdist

from functions.spectral import fixed_support_loss
from functions.storage import digest, read_json, save_json, verify_files


# 一、健康基底与损伤贡献；所有数组沿用冻结脑区顺序。
def loss_contributions(reference: dict, lost_counts: np.ndarray,
                       coverage: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """返回 135×129 的 G/W 节点贡献；排除节点为 NaN。"""
    coverage = np.asarray(coverage, dtype=float)
    if coverage.shape != (len(reference["active"]),) or not np.all(np.isfinite(coverage)):
        raise ValueError("灰质覆盖率必须为完整脑区的一维有限数组")
    if np.any((coverage < 0) | (coverage > 1)):
        raise ValueError("灰质覆盖率必须位于 [0,1]")
    local = coverage[:, None] * reference["vectors"] ** 2
    local[~reference["active"]] = np.nan
    _, _, network = fixed_support_loss(lost_counts, reference, contributions=True)
    return local, network


def combination_maps(reference: dict, lost_counts: np.ndarray, coverage: np.ndarray,
                     coefficients: np.ndarray, means: np.ndarray,
                     scales: np.ndarray) -> dict:
    """将标准化后的 258 个 G/W 系数变成两张脑图。

    顺序固定为 G1..G129,W1..W129。图为损伤贡献，不是带符号 u 的线性组合。
    两图节点求和 = 标准化影像得分 + 中心化常数。
    """
    modes = len(reference["values"])
    coefficients, means, scales = [np.asarray(value, dtype=float).reshape(-1)
                                    for value in (coefficients, means, scales)]
    if any(value.shape != (2 * modes,) for value in (coefficients, means, scales)):
        raise ValueError("组合系数、均值、标准差必须按 G 后 W 排成 258 维")
    if not all(np.isfinite(value).all() for value in (coefficients, means, scales)):
        raise ValueError("组合定义不能包含非有限值")
    if np.any(scales <= 0):
        raise ValueError("常量特征须使用标准差 1 和系数 0；不允许零标准差")
    local, network = loss_contributions(reference, lost_counts, coverage)
    features = np.r_[np.nansum(local, axis=0), np.nansum(network, axis=0)]
    raw = coefficients / scales  # 原始损伤量纲下的组合系数。
    local_map = local @ raw[:modes]
    network_map = network @ raw[modes:]
    offset = float(raw @ means)  # 从原始贡献总和转换到标准化得分时减去。
    score = float(coefficients @ ((features - means) / scales))
    error = float(np.nansum(local_map) + np.nansum(network_map) - offset - score)
    if not np.isclose(error, 0, atol=1e-9, rtol=1e-8):
        raise ValueError("空间贡献与固定组合得分不守恒")
    return {"g": local_map, "w": network_map, "score": score,
            "centering_constant": offset, "conservation_error": error}


def mean_combination_maps(reference: dict, files: Iterable[Path],
                          coefficients: np.ndarray, means: np.ndarray,
                          scales: np.ndarray) -> dict:
    """逐例读取紧凑 NPZ，在线平均贡献；不保存患者×脑区×模态大数组。"""
    local_sum = np.zeros(len(reference["active"]))
    network_sum = local_sum.copy()
    score_sum = 0.0
    count = 0
    offset = 0.0
    for path in files:
        with np.load(path) as patient:
            item = combination_maps(reference, patient["lost_counts"],
                                    patient["roi_coverage"], coefficients, means, scales)
        local_sum += np.nan_to_num(item["g"])
        network_sum += np.nan_to_num(item["w"])
        score_sum += item["score"]
        offset = item["centering_constant"]
        count += 1
    if not count:
        raise ValueError("组合空间图需要至少一个有效患者")
    local_sum /= count
    network_sum /= count
    local_sum[~reference["active"]] = np.nan
    network_sum[~reference["active"]] = np.nan
    return {"g": local_sum, "w": network_sum, "n": count,
            "score": score_sum / count, "centering_constant": offset,
            "object": "指定队列平均损伤贡献；受该队列病灶分布影响"}


def export_node_map(reference: dict, values: np.ndarray, directory: Path,
                    title: str, nifti: bool = False, figure: dict | None = None) -> None:
    """脑区表和展示 PNG；NIfTI 仅按需导出以降低存储。"""
    from functions.plotting import mode_png
    values = np.asarray(values, dtype=float).copy()
    if values.shape != reference["active"].shape:
        raise ValueError("空间图脑区数量不同")
    values[~reference["active"]] = np.nan
    directory.mkdir(parents=True, exist_ok=True)
    table = reference["nodes"][["index", "label", "group"]].copy()
    table["value"] = values
    table.to_csv(directory / "regions.tsv", sep="\t", index=False)
    figure = figure or {}
    mode_png(reference["template"], reference["atlas"], values, title,
             directory / "map.png", note="parcel loadings/contribution; excluded nodes masked",
             dpi=int(figure.get("dpi", 180)),
             cmap_signed=figure.get("colormap_signed", "coolwarm"),
             cmap_loss=figure.get("colormap_loss", "magma"))
    if nifti:
        import nibabel as nib
        mapped = np.r_[np.nan, values][reference["atlas"]].astype(np.float32)
        nib.save(nib.Nifti1Image(mapped, reference["template"].affine),
                 directory / "map.nii.gz")


def export_healthy_modes(reference: dict, directory: Path,
                         nifti: bool = False, figure: dict | None = None) -> dict:
    """全部非零健康模态的有符号背景图及平方载荷图；不是患者结果。"""
    directory.mkdir(parents=True, exist_ok=False)
    for position, mode in enumerate(reference["mode_ids"]):
        for kind in ("signed", "squared"):
            values = reference["vectors"][:, position].copy()
            if kind == "squared":
                values **= 2
            export_node_map(reference, values, directory / f"mode_{mode}" / kind,
                            f"Healthy mode {mode}: {kind}", nifti=nifti, figure=figure)
    record = {"modes": len(reference["mode_ids"]), "objects_per_mode": 2,
              "excluded_nodes": "NaN", "nifti": nifti,
              "mode_ids": [int(value) for value in reference["mode_ids"]],
              "figure": figure or {"dpi": 180, "colormap_signed": "coolwarm", "colormap_loss": "magma"}}
    save_json(directory / "manifest.json", record)
    return record


# 二、Moran 谱随机化。按半球保留均值、方差与给定空间权重下的 Moran I。
def spatial_weights(coordinates: np.ndarray, hemisphere: np.ndarray) -> np.ndarray:
    """体积脑区质心距离的倒数；不宣称这是表面测地距离或球面 spin。"""
    coordinates = np.asarray(coordinates, dtype=float)
    hemisphere = np.asarray(hemisphere)
    if coordinates.shape != (len(hemisphere), 3) or not np.isfinite(coordinates).all():
        raise ValueError("坐标必须是每个皮层脑区的有限三维 MNI 坐标")
    distance = cdist(coordinates, coordinates)
    same = hemisphere[:, None] == hemisphere[None, :]
    valid = same & ~np.eye(len(hemisphere), dtype=bool)
    if np.any(distance[valid] <= 0):
        raise ValueError("不同脑区不能共用完全相同的质心")
    return np.divide(1.0, distance, out=np.zeros_like(distance), where=valid)


def moran_surrogates(values: np.ndarray, weights: np.ndarray, hemisphere: np.ndarray,
                     n_perm: int = 10000, seed: int = 2026) -> tuple[np.ndarray, dict]:
    """singleton Moran 随机化；与 BrainSpace 的 MEM 符号翻转定义等价。

    使用与常数正交的完整子空间计算双中心化空间矩阵，避免数值零特征值误删。
    不做秩重排，因为秩重排会改变所保留的空间自相关。
    """
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    hemisphere = np.asarray(hemisphere)
    if values.ndim != 1 or not np.isfinite(values).all() or n_perm < 1:
        raise ValueError("零模型需要完整有限的一维皮层图和正整数重复次数")
    if weights.shape != (len(values), len(values)) or not np.allclose(weights, weights.T):
        raise ValueError("空间权重必须与脑图匹配且对称")
    random = np.random.default_rng(seed)
    nulls = np.empty((len(values), n_perm))
    diagnostics = []
    for label in np.unique(hemisphere):
        selected = np.flatnonzero(hemisphere == label)
        block = weights[np.ix_(selected, selected)]
        data = values[selected]
        if len(selected) < 3:
            raise ValueError("半球内有效脑区少于三个")
        centered = data - data.mean()
        basis = null_space(np.ones((1, len(selected))))
        _, vectors = np.linalg.eigh(basis.T @ block @ basis)
        mem = basis @ vectors
        coefficients = mem.T @ centered
        signs = random.choice(np.array([-1.0, 1.0]), size=(len(coefficients), n_perm))
        surrogate = mem @ (coefficients[:, None] * signs)
        nulls[selected] = surrogate + data.mean()
        variance_error = float(np.max(np.abs(np.sum(surrogate ** 2, axis=0)
                                             - centered @ centered)))
        spatial_error = float(np.max(np.abs(np.sum(surrogate * (block @ surrogate), axis=0)
                                            - centered @ block @ centered)))
        tolerance = 1e-8 * max(1.0, float(centered @ centered))
        if max(variance_error, spatial_error) > tolerance:
            raise ValueError("Moran 零模型没有保持方差/空间二次型")
        diagnostics.append({"hemisphere": str(label), "nodes": len(selected),
                            "variance_error": variance_error,
                            "spatial_quadratic_error": spatial_error})
    return nulls, {"method": "hemisphere_singleton_moran",
                   "distance": "inverse_euclidean_parcel_centroids_mm",
                   "n_perm": n_perm, "seed": seed, "hemispheres": diagnostics,
                   "limitation": "保留指定距离权重的 Moran I，未保留完整经验分布；不是球面旋转"}


def correlation_test(values: np.ndarray, annotations: np.ndarray,
                     nulls: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """双侧空间 p=(1+极端零模型数量)/(1+重复次数)，不使用独立脑区 t 检验。"""
    values = np.asarray(values, dtype=float)
    annotations = np.asarray(annotations, dtype=float)
    centered = values - values.mean()
    annotated = annotations - annotations.mean(axis=0)
    denominator = np.linalg.norm(centered) * np.linalg.norm(annotated, axis=0)
    if np.any(denominator <= 1e-15):
        raise ValueError("常量脑图或常量注释无法计算相关")
    observed = np.clip(centered @ annotated / denominator, -1, 1)
    random = nulls - nulls.mean(axis=0)
    simulated = random.T @ annotated / (
        np.linalg.norm(random, axis=0)[:, None] * np.linalg.norm(annotated, axis=0)[None, :])
    probability = (1 + (np.abs(simulated) >= np.abs(observed)[None, :] - 1e-12).sum(axis=0)) / (1 + nulls.shape[1])
    return observed, probability


def family_bh(probabilities: np.ndarray) -> np.ndarray:
    """失败检验保留计划检验数；NA 当作 p=1 排序，最终仍输出 NA。"""
    values = np.asarray(probabilities, dtype=float)
    valid = np.isfinite(values)
    filled = np.where(valid, values, 1.0)
    order = np.argsort(filled, kind="stable")
    ordered = filled[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(ordered[::-1])[::-1].clip(0, 1)
    result = np.empty(len(values))
    result[order] = adjusted
    result[~valid] = np.nan
    return result


# 三、正式资源验签与阶段入口。
def load_annotations(root: Path, resource_path: str = "reference/annotation") -> dict:
    directory = (root / resource_path).resolve()
    directory.relative_to(root.resolve())
    manifest = read_json(directory / "manifest.json")
    verify_files(directory, manifest["files"])
    if digest(root / "reference/nodes.tsv") != manifest["project_nodes_sha256"]:
        raise ValueError("脑区表发生改变，需要重新核对注释空间")
    if digest(root / "reference/atlas.nii.gz") != manifest["project_atlas_sha256"]:
        raise ValueError("脑区图谱发生改变，需要重新核对注释空间")
    if manifest["alignment"]["status"] != "passed":
        raise ValueError("注释分区对应尚未通过")
    nodes = pd.read_csv(directory / "cortical_nodes.tsv", sep="\t")
    result = {"manifest": manifest, "nodes": nodes,
              "hansen": pd.read_csv(directory / "hansen.tsv", sep="\t"),
              "neurosynth": pd.read_csv(directory / "neurosynth.tsv", sep="\t")}
    if nodes["index"].to_list() != list(range(1, 101)):
        raise ValueError("冻结皮层节点编号必须为1至100且顺序一致")
    for family in ("hansen", "neurosynth"):
        if result[family]["index"].to_list() != nodes["index"].to_list():
            raise ValueError(f"{family} 脑区编号与核对表不同")
    return result


def run_annotation(root: Path, config: dict, input_dir: Path, output_dir: Path) -> dict:
    """输入 maps_manifest.json 的预设 G/W 图，输出两个独立家族的空间注释。

    maps 条目需要 pattern_id、channel、map_file、validated。
    map_file 为相对 input_dir 路径，TSV 必须包含 index 与 value，完整135行。
    未验证对象保留 NA 结果及原因，不产生机制结论。
    """
    options = config.get("annotation", config.get("annotations", {}))
    resources = load_annotations(root, options.get("resources", "reference/annotation"))
    if options.get("null_method", "hemisphere_singleton_moran") != "hemisphere_singleton_moran":
        raise ValueError("本模块只实现 hemisphere_singleton_moran，不能忽略未实现的零模型配置")
    n_perm = int(options.get("permutations", options.get("n_perm", 10000)))
    seed = int(options.get("seed", 2026))
    fdr = options.get("fdr", {})
    if fdr.get("method", "BH") != "BH":
        raise ValueError("正式空间注释只实现预设的BH校正")
    q_threshold = float(fdr.get("q", 0.05))
    if not 0 < q_threshold < 1:
        raise ValueError("FDR阈值应位于0和1之间")
    planned = read_json(input_dir / "maps_manifest.json")["maps"]
    if len({(item["pattern_id"], item["channel"]) for item in planned}) != len(planned):
        raise ValueError("注释对象清单存在重复组合/通道")
    output_dir.mkdir(parents=True, exist_ok=True)
    nodes = resources["nodes"]
    weights = spatial_weights(nodes[["x_mm", "y_mm", "z_mm"]].to_numpy(), nodes.hemisphere.to_numpy())
    records, diagnostics = [], []
    for item in planned:
        if item["channel"] not in {"g", "w"}:
            raise ValueError("正式注释只接受 g、w 两个通道")
        reason, values, nulls = "", None, None
        if not item.get("validated", False):
            reason = "validation_not_passed"
        else:
            path = (input_dir / item["map_file"]).resolve()
            path.relative_to(input_dir.resolve())
            table = pd.read_csv(path, sep="\t")
            if len(table) != 135 or table["index"].duplicated().any() or set(table["index"]) != set(range(1, 136)):
                raise ValueError("注释图必须给出完整、唯一的135个脑区编号")
            values = table.set_index("index").loc[nodes["index"], "value"].to_numpy()
            if not np.isfinite(values).all():
                raise ValueError("主分析100个皮层节点必须全部有有效覆盖，不能静默删除")
            if np.std(values) <= 1e-15:
                reason = "constant_map"
            else:
                nulls, diagnostic = moran_surrogates(values, weights, nodes.hemisphere.to_numpy(), n_perm, seed)
                diagnostic.update({"pattern_id": item["pattern_id"], "channel": item["channel"]})
                diagnostics.append(diagnostic)
        for family in ("neurosynth", "hansen"):
            table = resources[family]
            names = list(table.columns[1:])
            correlations, probabilities = (np.full(len(names), np.nan), np.full(len(names), np.nan))
            available = np.isfinite(table.iloc[:, 1:].to_numpy()).all(axis=0)
            if not reason:
                correlations[available], probabilities[available] = correlation_test(values, table.iloc[:, 1:].to_numpy()[:, available], nulls)
            for position, (name, correlation, probability) in enumerate(zip(names, correlations, probabilities)):
                records.append({"pattern_id": item["pattern_id"], "channel": item["channel"],
                                "family": family, "annotation": name, "r": correlation,
                                "p": probability, "status": reason or ("estimated" if available[position] else "resource_unavailable"), "nodes": 100})
    result = pd.DataFrame(records)
    if len(result):
        result["q"] = np.nan
        for family in ("neurosynth", "hansen"):
            selected = result.family == family
            result.loc[selected, "q"] = family_bh(result.loc[selected, "p"].to_numpy())
        result["significant"] = result["q"].le(q_threshold)
        result.to_csv(output_dir / "associations.tsv", sep="\t", index=False)
        annotation_png(result, output_dir / "associations.png", int(options.get("figure", {}).get("dpi", 180)))
    save_json(output_dir / "null_diagnostics.json", {"maps": diagnostics})
    summary = {"status": "completed" if len(result) else "no_eligible_maps",
               "planned_maps": len(planned), "estimated_maps": len(diagnostics),
               "tests": len(result), "n_perm": n_perm, "scope": "100 cortical parcels only",
               "fdr_q_threshold": q_threshold,
               "other_structures": "皮层下、小脑、脑干不进入首版空间统计",
               "resource_status": resources["manifest"]["status"],
               "unavailable_hansen": sorted(set(resources["manifest"]["hansen_planned"])
                                            - set(resources["manifest"]["hansen_available"])),
               "resource_id": resources["manifest"]["resource_id"]}
    save_json(output_dir / "summary.json", summary)
    return summary


def annotation_png(table: pd.DataFrame, path: Path, dpi: int = 180) -> None:
    """完整检验表保留所有术语；图只展示各资源绝对相关最大的20项。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(14, 8), dpi=dpi)
    for axis, family in zip(axes, ("neurosynth", "hansen")):
        selected = table.loc[(table.family == family) & table.r.notna()].copy()
        selected["absolute_r"] = selected.r.abs()
        selected = selected.nlargest(20, "absolute_r").sort_values("r")
        labels = selected.pattern_id + ":" + selected.channel + " | " + selected.annotation
        axis.barh(range(len(selected)), selected.r, color=np.where(selected.r >= 0, "#a94442", "#376c9b"))
        axis.set_yticks(range(len(selected)), labels, fontsize=7)
        axis.set(xlim=(-1, 1), xlabel="Spatial Pearson r", title=family)
        axis.axvline(0, color="black", linewidth=.5)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)
