"""配置读取与显式数据合同；不解释任何 Python/R 表达式。"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import yaml


def validate_contract(config: dict) -> None:
    """固定方法属性不能伪装成开关；改变这些属性须修改方法并重新验证代码。"""
    fixed = {
        "paths.reference": "reference",
        "imaging.space": "MNI152NLin6Asym",
        "imaging.voxel_size_mm": 1.0,
        "imaging.healthy_contract": {"nodes":135,"active_nodes":130,"nonzero_modes":129,"edges":1375},
        "clinical.merge.unmatched_policy": "report",
        "clinical.merge.conflict_policy": "error",
        "clinical.death.policy": "mrs_equals_6",
        "statistics.channels": ["g", "w"],
        "statistics.fdr.method": "BH",
        "statistics.fdr.failed_tests_in_denominator": True,
        "statistics.mi.outcome_policy": "observed_only",
        "statistics.mi.inference_outcome_predictor": True,
        "statistics.mi.learning_outcome_predictor": False,
        "statistics.mi.save_completed_images": False,
        "statistics.models.ordinal.package": "ordinal",
        "statistics.models.ordinal.link": "logit",
        "statistics.models.linear.covariance": "HC3",
        "statistics.models.longitudinal.time": "categorical",
        "statistics.models.longitudinal.random_intercept": "subject_id",
        "statistics.models.joint_test": "D1",
        "statistics.models.coefficient_pooling": "Rubin",
        "patterns.package": "ordinalNet",
        "patterns.standardize_within_training": True,
        "patterns.imputation_uses_outcome": False,
        "patterns.covariates_penalized": False,
        "patterns.score_contains_covariates": False,
        "validation.forbid_subject_overlap": True,
        "validation.require_replication_for_annotation": True,
        "validation.inference_mi_uses_observed_outcome": True,
        "annotation.cortex_only": True,
        "annotation.extend_subcortex": False,
        "annotation.figure.format": "png",
        "annotation.figure.dpi": 300,
        "annotation.fdr.method": "BH",
        "annotation.fdr.families": ["neurosynth", "hansen"],
        "runtime.checkpoint": True,
        "runtime.atomic_commit": True,
        "runtime.compact_outputs": True,
        "runtime.patient_qc_images": False,
        "runtime.manual_qc_tables": False,
    }
    for dotted, expected in fixed.items():
        actual = config
        for key in dotted.split("."):
            if key not in actual:
                raise ValueError(f"缺少固定方法属性 {dotted}")
            actual = actual[key]
        if actual != expected:
            raise ValueError(f"{dotted} 为固定方法属性，必须是 {expected!r}，不能静默忽略修改")
    mi = config["statistics"]["mi"]
    if int(mi["m"]) < 2 or int(mi["maxit"]) < 2 or int(mi["increment"]) < 1 or int(mi["max_iterations"]) < int(mi["maxit"]):
        raise ValueError("MI 至少两套，迭代至少两次、递增为正且上限不小于初始值")
    patterns = config["patterns"]
    if int(patterns["folds"]) < 2 or int(patterns["n_lambda"]) < 2:
        raise ValueError("组合学习至少两折且至少两个惩罚参数")
    if not 0 <= float(patterns["alpha"]) <= 1 or not 0 < float(patterns["lambda_min"]) < float(patterns["lambda_max"]):
        raise ValueError("alpha 必须在 [0,1]，且 0 < lambda_min < lambda_max")
    if patterns["selection"] not in {"one_standard_error", "one_se"}:
        raise ValueError("预设组合选择规则为一标准误规则")
    if not 0 < float(config["statistics"]["fdr"]["q"]) < 1:
        raise ValueError("FDR 阈值必须在 (0,1)")
    allowed = {"preflight","lesion","clinical","imputation","association","patterns","validation","annotation","report"}
    stages = config["runtime"]["stages"]
    if not stages or len(stages) != len(set(stages)) or set(stages)-allowed:
        raise ValueError("runtime.stages 包含未知或重复阶段")


def load_config(path: Path | str) -> dict:
    """读取主配置，核查结构，并为旧病灶总控保留只读兼容别名。"""
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError("配置必须为键值结构")
    if "imaging" not in config:  # 仅供历史三例等价测试使用。
        return config
    required = ["paths", "subjects", "clinical", "statistics", "patterns", "validation", "annotation", "runtime"]
    missing = [name for name in required if not isinstance(config.get(name), dict)]
    if missing:
        raise ValueError(f"缺少配置分区：{missing}")
    # 发现集与验证集分别检查文件模式和ID规则；正式配置不能退回单一病灶目录。
    lesion_paths = config["paths"].get("lesion")
    if not isinstance(lesion_paths, dict) or set(lesion_paths) != {"discovery", "validation"}:
        raise ValueError("paths.lesion 必须同时指定 discovery 和 validation")
    clinical_paths = config["paths"].get("clinical")
    if not isinstance(clinical_paths, dict) or set(clinical_paths) != {"discovery", "validation"}:
        raise ValueError("paths.clinical 必须同时指定 discovery 和 validation")
    cohort_settings = config["subjects"].get("cohorts")
    if not isinstance(cohort_settings, dict) or set(cohort_settings) != {"discovery", "validation"}:
        raise ValueError("subjects.cohorts 必须同时指定 discovery 和 validation")
    for cohort_name in ("discovery", "validation"):
        settings = cohort_settings[cohort_name]
        pattern = str(settings.get("pattern", ""))
        id_rule = settings.get("id_rule", {}) or {}
        id_regex = settings.get("id_regex") or id_rule.get("id_regex")
        if id_regex:
            import re
            try:
                compiled = re.compile(str(id_regex))
            except re.error as error:
                raise ValueError(f"{cohort_name} id_rule.id_regex 无效: {error}") from error
            if "subject_id" not in compiled.groupindex:
                raise ValueError(f"{cohort_name} id_rule.id_regex 必须包含命名组 subject_id")
        elif pattern.count("*") != 1 or any(value in pattern for value in ["?", "[", "]", "/", "\\"]):
            raise ValueError(f"subjects.cohorts.{cohort_name}.pattern 必须恰好包含一个 *，不能含其他通配符或路径")
        elif not id_rule.get("prefix") and not id_rule.get("suffix"):
            raise ValueError(f"{cohort_name} 必须提供 id_rule.prefix/suffix 或 id_regex")
    runtime = config["runtime"]
    if not 1 <= int(runtime["workers"]) <= 32 or int(runtime["threads"]) < 1:
        raise ValueError("workers 必须为 1–32，threads 必须为正整数")
    if set(config["statistics"]["families"]) != {"nihss", "mrs_single", "mrs_longitudinal", "moca"}:
        raise ValueError("正式统计必须显式保留四个预设检验家族（MoCA 3/6/12月共用一个家族）")
    variables = config["clinical"]["variables"]
    allowed = {"continuous", "binary", "nominal", "ordinal"}
    for name, specification in variables.items():
        if specification.get("type") not in allowed:
            raise ValueError(f"变量 {name} 的类型无效")
        if specification["type"] in {"binary", "nominal", "ordinal"} and not specification.get("levels"):
            raise ValueError(f"分类变量 {name} 必须明确 levels")
    if len(config["patterns"]["targets"]) != 20 or len(set(config["patterns"]["targets"])) != 20:
        raise ValueError("组合学习固定为 20 个不同目标")
    validate_contract(config)
    config = deepcopy(config)
    config["space"] = config["imaging"]["space"]
    config["maximum_volume_change"] = config["imaging"]["maximum_volume_change"]
    config["workers"] = runtime["workers"]
    config["threads"] = runtime["threads"]
    config["save_jmax"] = config["imaging"].get("save_jmax", True)
    config["participants"] = config["subjects"].get("participants")
    config["cohorts"] = config["subjects"]["cohorts"]
    return config
