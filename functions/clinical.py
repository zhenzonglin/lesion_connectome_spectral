"""临床数据合同：SAS 导入、明确 ID 合并、编码核查、派生变量与影像连接。

本模块不插补、不拟合模型。所有原始输入只读，核查表写入调用者指定输出目录。
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import re
import numpy as np
import pandas as pd
from functions.storage import digest, save_json


# 一、配置与输入核查
def clinical_directory(root: Path, config: dict, cohort: str | None = None) -> Path:
    """返回队列临床目录；支持 paths.clinical 的双队列字典。"""
    configured = config["paths"]["clinical"]
    if isinstance(configured, dict):
        if cohort is None:
            raise ValueError("双队列临床路径需要指定 cohort")
        configured = configured.get(cohort)
        if not configured:
            raise ValueError(f"paths.clinical 缺少 {cohort} 路径")
    elif cohort is not None and config.get("clinical", {}).get("cohorts", {}).get(cohort) is not None:
        configured = str(configured) + f"/{cohort}"
    directory = (root / configured).resolve()
    directory.relative_to(root.resolve())
    return directory


def source_path(root: Path, config: dict, filename: str, cohort: str | None = None) -> Path:
    """SAS 和映射表只允许位于对应 clinical 目录内。"""
    directory = clinical_directory(root, config, cohort)
    path = (directory / filename).resolve()
    path.relative_to(directory)
    return path


def _sas_files(root: Path, config: dict, cohort: str) -> list[Path]:
    directory = clinical_directory(root, config, cohort)
    # SAS 扩展名大小写不影响发现；目录中的其他文件不参与读取。
    return sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() == ".sas7bdat"
    )


def clinical_input_files(root: Path, config: dict) -> list[Path]:
    """自动识别 discovery/validation 子目录中的全部 SAS 文件。"""
    files = []
    for cohort in ["discovery", "validation"]:
        try:
            files.extend(_sas_files(root, config, cohort))
        except (ValueError, OSError):
            continue
        settings = config.get("clinical", {}).get("cohorts", {}).get(cohort, {})
        for key in ["member_file", "mapping_file"]:
            filename = settings.get(key)
            if filename:
                files.append(source_path(root, config, filename, cohort))
        for source in settings.get("sources", []):
            if source.get("file"):
                files.append(source_path(root, config, source["file"], cohort))
            if source.get("id", {}).get("mapping_file"):
                files.append(source_path(root, config, source["id"]["mapping_file"], cohort))
    return sorted(set(files))


def _id_fields(settings: dict) -> list[str]:
    identity = settings.get("id", {}) or {}
    fields = settings.get("id_fields") or identity.get("fields")
    if fields:
        return [str(value) for value in fields]
    field = settings.get("id_field") or settings.get("id_column") or identity.get("field") or identity.get("column")
    return [str(field)] if field else []


def _field_for_cohort(config: dict, cohort: str, name: str) -> dict:
    specification = deepcopy(config["clinical"]["variables"].get(name, {}))
    settings = config["clinical"]["cohorts"].get(cohort, {})
    overrides = settings.get("variable_overrides", {}).get(name)
    if overrides is None:
        overrides = settings.get("field_map", {}).get(name)
    if isinstance(overrides, str):
        overrides = {"field": overrides}
    if overrides:
        specification.update(overrides)
    return specification


def outcome_definitions(config: dict) -> list[dict]:
    """生成临床结局定义；家族由 config.statistics.families 固定。"""
    definitions = []
    variables = config["clinical"]["variables"]
    for family, specification in config["statistics"]["families"].items():
        if not specification.get("enabled", True):
            continue
        for name in specification.get("outcomes", []):
            variable = variables.get(name, {"type": "ordinal", "levels": list(range(7)), "title": name})
            definitions.append({
                "name": name, "family": family, "model": specification.get("model"),
                **{key: value for key, value in variable.items()
                   if key in ["type", "levels", "title"]},
                "covariates": config["clinical"]["covariate_sets"][specification["covariate_set"]],
            })
    return definitions


def required_variables(config: dict) -> set[str]:
    """追溯启用结局、协变量、组合目标及派生变量的输入。"""
    required = set()
    for outcome in outcome_definitions(config):
        if outcome["name"] == "mrs_longitudinal":
            required.update(config["statistics"]["families"]["mrs_longitudinal"].get("visits", {}).keys())
        else:
            required.add(outcome["name"])
        required.update(outcome["covariates"])
    if config.get("patterns", {}).get("enabled", True):
        required.update(config["patterns"].get("targets", []))
        for name in config["patterns"].get("targets", []):
            covariate_key = "mrs" if name.startswith("mrs_") else "nihss"
            if covariate_key in config["clinical"]["covariate_sets"]:
                required.update(config["clinical"]["covariate_sets"][covariate_key])
    for _ in range(len(config["clinical"].get("derived", {})) + 1):
        previous = required.copy()
        for name, specification in config["clinical"].get("derived", {}).items():
            if name in required:
                required.update(specification.get("inputs", []))
        if previous == required:
            break
    return required


def _required_fields(config: dict, cohort: str) -> set[str]:
    required = required_variables(config)
    fields = set()
    for name in required:
        specification = _field_for_cohort(config, cohort, name)
        if specification.get("source") not in {"derived", "imaging"} and specification.get("field"):
            fields.add(str(specification["field"]))
    return fields


def preflight_clinical(root: Path, config: dict, cohort: str = "discovery") -> list[str]:
    """自动扫描 SAS 并检查 ID；discovery 字段严格，validation 允许缺少非共同结局。"""
    issues = []
    settings = config.get("clinical", {}).get("cohorts", {}).get(cohort)
    if settings is None:
        return [f"未配置队列 {cohort}"]
    # 临床字段合同需要人工确认后才能进入正式阶段。
    if not settings.get("confirmed", False):
        issues.append(f"{cohort}: clinical.cohorts.{cohort}.confirmed 必须为 true")
    rules_confirmed = settings.get("rules_confirmed", {}) or {}
    if not rules_confirmed.get("id_mapping", False):
        issues.append(f"{cohort}: clinical.cohorts.{cohort}.rules_confirmed.id_mapping 必须为 true")
    if not rules_confirmed.get("mrs_death", False):
        issues.append(f"{cohort}: clinical.cohorts.{cohort}.rules_confirmed.mrs_death 必须为 true")
    try:
        files = _sas_files(root, config, cohort)
    except (KeyError, ValueError, OSError) as error:
        return [f"{cohort}: 临床目录无效：{error}"]
    if not files:
        issues.append(f"{cohort}: clinical 目录没有 .sas7bdat 文件")
        return issues
    fields = _id_fields(settings)
    if not fields:
        issues.append(f"{cohort}: 必须指定 id_field 或 id_fields")
    # 读取元数据列名，避免等待完整临床阶段才发现 ID 缺失。
    available = set()
    for path in files:
        try:
            frame, _ = read_sas(path)
            available.update(map(str, frame.columns))
        except Exception as error:
            issues.append(f"{cohort}: 无法读取 {path.name}：{error}")
    missing_ids = set(fields) - available
    if missing_ids:
        issues.append(f"{cohort}: SAS 缺少 ID 字段：{sorted(missing_ids)}")
    missing_fields = _required_fields(config, cohort) - available
    if cohort == "discovery" and missing_fields:
        issues.append(f"{cohort}: SAS 缺少发现集必需字段：{sorted(missing_fields)}")
    settings.setdefault("rules_confirmed", {})
    return issues


# 二、ID 与 SAS 导入：绝不默认取第一条重复记录
def normalize_ids(values: pd.Series, specification: dict, mapping: pd.DataFrame | None = None) -> pd.Series:
    """将显式文本/数值规则转为稳定 ID；不隐式去 ses、改大小写或补零。"""
    if values.isna().any():
        raise ValueError("ID 存在缺失")
    if specification["kind"] == "numeric":
        numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
        if not np.isfinite(numeric).all() or not np.equal(numeric, np.floor(numeric)).all() or (np.abs(numeric) > 2**53).any():
            raise ValueError("数值 ID 必须是可精确表示的整数")
        identifiers = pd.Series([str(int(value)) for value in numeric], index=values.index)
        if mapping is None:
            width = specification.get("width")
            if not isinstance(width, int) or width < 1:
                raise ValueError("数值 ID 必须明确正整数补零宽度")
            identifiers = identifiers.str.zfill(width)
    else:
        if any(not isinstance(value, str) for value in values):
            raise ValueError("id.kind=string 但 SAS ID 不是字符串；需明确配置 numeric")
        identifiers = values.copy()
        if specification.get("strip", False):
            identifiers = identifiers.str.strip()
    if mapping is not None:
        if set(mapping.columns) != {"source_id", "subject_id"} or mapping["source_id"].duplicated().any() or mapping["subject_id"].duplicated().any():
            raise ValueError("ID 对应表必须有 source_id、subject_id 两列且均唯一")
        identifiers = identifiers.map(mapping.set_index("source_id")["subject_id"])
        if identifiers.isna().any():
            raise ValueError("ID 对应表缺少 SAS 编号")
    if not identifiers.str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*").all():
        raise ValueError("ID 含空白、空字符串或不安全字符；请显式修正 ID 规则")
    return identifiers.astype(str)


def read_sas(path: Path) -> tuple[pd.DataFrame, dict]:
    """保留 SAS 特殊缺失元信息后显式转缺失；不使用显示标签替代数值编码。"""
    import pyreadstat
    frame, metadata = pyreadstat.read_sas7bdat(str(path), user_missing=True, disable_datetime_conversion=True)
    special = getattr(metadata, "missing_user_values", {}) or {}
    for column, missing in special.items():
        if column in frame:
            frame[column] = frame[column].mask(frame[column].isin(missing))
    information = {"column_labels": getattr(metadata, "column_names_to_labels", {}), "special_missing": special,
                   "file_label": getattr(metadata, "file_label", None), "number_rows": len(frame)}
    return frame, information


def write_audit(output_dir: Path, name: str, records: list[dict], columns: list[str] | None = None) -> None:
    directory = output_dir / "audit"
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=columns).to_csv(directory / f"{name}.tsv", sep="\t", index=False)


def _normalize_composite_ids(frame: pd.DataFrame, settings: dict) -> pd.Series:
    """按单列或多列拼接生成稳定 subject_id；数值字段必须显式配置格式。"""
    fields = _id_fields(settings)
    if not fields or any(field not in frame for field in fields):
        raise ValueError(f"临床表缺少 ID 字段：{fields}")
    identity = settings.get("id", {}) or {}
    if len(fields) == 1:
        specification = dict(identity)
        specification.setdefault("kind", settings.get("id_kind", "string"))
        specification.setdefault("strip", True)
        return normalize_ids(frame[fields[0]], specification)
    field_specs = identity.get("field_specs", {}) or {}
    shared_kind = identity.get("kind", "string")
    shared_width = identity.get("width")
    normalized_parts = []
    for field in fields:
        specification = dict(field_specs.get(field, {}))
        if not specification:
            specification = {"kind": shared_kind}
            if shared_width is not None:
                specification["width"] = shared_width
            specification["strip"] = True
        elif "kind" not in specification:
            specification["kind"] = shared_kind
        # 不能把 SAS 数值列静默转成 12.0；需设置 field_specs[field].kind/width。
        normalized_parts.append(normalize_ids(frame[field], specification))

    separator = str(settings.get("id_separator", identity.get("separator", "")))
    result = normalized_parts[0].astype(str)
    for part in normalized_parts[1:]:
        result = result + separator + part.astype(str)
    if not result.str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*").all():
        raise ValueError("复合 ID 含不安全字符")
    return result


def load_sources(root: Path, config: dict, cohort: str, output_dir: Path) -> tuple[pd.DataFrame, list[dict]]:
    """自动读取队列目录的全部 SAS；同名列仅在非冲突时合并。"""
    settings = config["clinical"]["cohorts"][cohort]
    files = _sas_files(root, config, cohort)
    if not files:
        raise ValueError(f"{cohort}: 没有可读取的 SAS 文件")
    merged = None
    provenance, conflicts = [], []
    for path in files:
        frame, metadata = read_sas(path)
        frame = frame.copy()
        frame["subject_id"] = _normalize_composite_ids(frame, settings)
        if frame.subject_id.duplicated().any():
            duplicates = frame.loc[frame.subject_id.duplicated(keep=False), "subject_id"].tolist()
            write_audit(output_dir, "duplicate_records",
                        [{"source": path.name, "subject_id": value} for value in duplicates])
            raise ValueError(f"SAS 文件 {path.name} 存在重复 subject_id")
        frame = frame.drop(columns=_id_fields(settings), errors="ignore")
        if merged is None:
            merged = frame
        else:
            left = merged.set_index("subject_id")
            right = frame.set_index("subject_id")
            overlap = set(left.index) & set(right.index)
            for column in set(left.columns) & set(right.columns):
                both = left.loc[left.index.intersection(right.index), column]
                other = right.loc[both.index, column]
                conflict = both.notna() & other.notna() & (both.astype(str) != other.astype(str))
                if conflict.any():
                    conflicts.extend({"source": path.name, "variable": column,
                                      "subject_id": str(index)}
                                     for index in both.index[conflict])
            if conflicts:
                write_audit(output_dir, "source_conflicts", conflicts,
                            ["source", "variable", "subject_id"])
                raise ValueError("多个 SAS 文件同名字段存在冲突，见 audit/source_conflicts.tsv")
            merged = left.combine_first(right).reset_index()
            unmatched = [{"source": path.name, "subject_id": value, "side": "previous_only"}
                         for value in sorted(set(left.index)-set(right.index))]
            unmatched += [{"source": path.name, "subject_id": value, "side": "new_only"}
                          for value in sorted(set(right.index)-set(left.index))]
            if unmatched:
                write_audit(output_dir, "source_unmatched", unmatched,
                            ["source", "subject_id", "side"])
        provenance.append({"name": path.stem, "file": str(path.relative_to(root)),
                           "sha256": digest(path), "operation": "auto_join", **metadata})
    if merged is None:
        raise ValueError("没有临床来源")
    return merged, provenance


# 三、字段映射、缺失与临床有效性
def select_field(frame: pd.DataFrame, specification: dict, name: str) -> pd.Series:
    field = specification.get("field")
    source = specification.get("source")
    if field in frame:
        return frame[field].copy()
    candidates = [f"{source}::{field}"] if source else [column for column in frame if column.endswith("::" + str(field))]
    candidates = [column for column in candidates if column in frame]
    if len(candidates) != 1:
        raise ValueError(f"变量 {name} 字段 {field} 有 {len(candidates)} 个来源；须明确 source，不能默认取第一列")
    return frame[candidates[0]].copy()


def map_variables(raw: pd.DataFrame, config: dict, cohort: str, output_dir: Path, allow_missing: bool = False) -> pd.DataFrame:
    """按字段独立编码，异常保留审计并报错；未启用字段可缺失。"""
    variables = {
        name: _field_for_cohort(config, cohort, name)
        for name in config["clinical"]["variables"]
    }
    required = required_variables(config)
    result = raw[["subject_id"]].copy()
    invalid = []
    for name, specification in variables.items():
        if specification.get("source") in {"derived", "imaging"}:
            continue
        if name not in required:
            result[name] = np.nan
            continue
        try:
            values = select_field(raw, specification, name)
        except ValueError as error:
            if name in required and not allow_missing:
                write_audit(output_dir, "field_conflicts", [{"variable": name, "reason": str(error)}])
                raise
            write_audit(output_dir, "missing_fields",
                        [{"variable": name, "reason": str(error), "cohort": cohort}])
            result[name] = np.nan
            continue
        values = values.mask(values.isin(specification.get("missing", [])))
        if specification.get("recode"):
            recoding = {str(key): value for key, value in specification["recode"].items()}
            def recode(value: object) -> object:
                if pd.isna(value):
                    return np.nan
                key = str(int(value)) if isinstance(value, (float, np.floating)) and value.is_integer() else str(value)
                return recoding.get(key, value)
            values = values.map(recode)
        numeric_levels = all(isinstance(value, (int, float)) for value in specification.get("levels", []))
        if specification["type"] == "continuous" or numeric_levels:
            numeric = pd.to_numeric(values, errors="coerce")
            bad = values.notna() & numeric.isna()
            bad |= numeric.notna() & ~np.isfinite(numeric)
        else:
            numeric = values.copy()
            bad = pd.Series(False, index=values.index)
        if "levels" in specification:
            bad |= numeric.notna() & ~numeric.isin(specification["levels"])
        if "range" in specification:
            lower, upper = specification["range"]
            bad |= numeric.notna() & ~numeric.between(lower, upper)
        invalid.extend({"subject_id": raw.loc[index, "subject_id"], "variable": name, "raw_value": str(values.loc[index]), "reason": "invalid_code_or_range"} for index in raw.index[bad])
        result[name] = numeric
    write_audit(output_dir, "invalid_values", invalid, ["subject_id", "variable", "raw_value", "reason"])
    if invalid:
        raise ValueError("临床编码或范围异常，见 audit/invalid_values.tsv；未静默转为缺失")
    return result


def derive_variables(frame: pd.DataFrame, derived: dict) -> pd.DataFrame:
    """允许的派生运算只有 log1p、完整组成项求和；禁用 eval。"""
    result = frame.copy()
    pending = deepcopy(derived)
    while pending:
        progressed = False
        for name, specification in list(pending.items()):
            inputs = specification["inputs"]
            if not all(column in result for column in inputs):
                continue
            if specification["operation"] == "sum_complete":
                result[name] = result[inputs].sum(axis=1, min_count=len(inputs))
            elif specification["operation"] == "log1p" and len(inputs) == 1:
                if (result[inputs[0]].dropna() < 0).any():
                    raise ValueError("log1p 体积不能为负")
                result[name] = np.log1p(result[inputs[0]])
            else:
                raise ValueError(f"不支持的派生运算：{specification['operation']}")
            del pending[name]
            progressed = True
        if not progressed:
            raise ValueError(f"派生变量缺少输入或存在循环：{list(pending)}")
    return result


def check_death(raw: pd.DataFrame, mapped: pd.DataFrame, config: dict, output_dir: Path, cohort: str = "discovery") -> pd.DataFrame:
    """仅由对应月份 mRS=6 推断死亡；死后 MoCA 作为结构性缺失。"""
    result = mapped.copy()
    conflicts = []
    required = required_variables(config)
    previous_dead = pd.Series(False, index=result.index)
    months = sorted(
        int(re.search(r"_(\d+)m$", name).group(1))
        for name in required
        if re.fullmatch(r"mrs_\d+m", name)
    )
    for month in months:
        mrs, moca = f"mrs_{month}m", f"moca_{month}m"
        if mrs not in result:
            continue
        values = pd.to_numeric(result[mrs], errors="coerce")
        dead_now = values.eq(6)
        alive_after_death = previous_dead & values.notna() & ~dead_now
        for index in result.index[alive_after_death]:
            conflicts.append({"subject_id": result.loc[index, "subject_id"],
                              "month": month, "reason": "mrs_observed_after_prior_death"})
        dead = previous_dead | dead_now
        result[f"death_{month}m"] = dead.astype(float)
        if moca in result:
            for index in result.index[dead & result[moca].notna()]:
                conflicts.append({"subject_id": result.loc[index, "subject_id"],
                                  "month": month, "reason": "moca_observed_after_death"})
            result.loc[dead, moca] = np.nan
            result[f"{moca}_structural_missing"] = dead
        previous_dead = dead
    write_audit(output_dir, "death_conflicts", conflicts, ["subject_id", "month", "reason"])
    if conflicts:
        raise ValueError("死亡状态与临床结局冲突，见 audit/death_conflicts.tsv")
    return result


def assert_disjoint(discovery_ids: set[str], validation_ids: set[str]) -> None:
    overlap = discovery_ids & validation_ids
    if overlap:
        raise ValueError(f"发现集与验证集患者重叠：{sorted(overlap)[:20]}（共{len(overlap)}例）")


# 四、影像与临床整合；输出固定宽表供 R 推断和组合学习复用
def prepare_clinical(root: Path, config: dict, feature_dir: Path, output_dir: Path, cohort: str = "discovery") -> dict:
    issues = preflight_clinical(root, config, cohort)
    output_dir.mkdir(parents=True, exist_ok=True)
    if issues:
        save_json(output_dir / "preflight.json", {"status": "blocked_missing_input", "issues": issues})
        raise ValueError("；".join(issues))
    raw, provenance = load_sources(root, config, cohort, output_dir)
    settings = config["clinical"]["cohorts"][cohort]
    if settings.get("member_file"):
        members = pd.read_csv(source_path(root, config, settings["member_file"], cohort), sep="\t", dtype=str)
        if list(members.columns) != ["subject_id"] or members.subject_id.duplicated().any() or members.subject_id.isna().any():
            raise ValueError("成员表必须只有唯一 subject_id 列")
        absent = set(members.subject_id) - set(raw.subject_id)
        if absent:
            raise ValueError(f"成员表 {len(absent)} 人没有临床来源")
        raw = raw[raw.subject_id.isin(members.subject_id)].reset_index(drop=True)
    mapped = check_death(raw, map_variables(raw, config, cohort, output_dir, allow_missing=(cohort == "validation")), config, output_dir, cohort)
    features = pd.read_csv(feature_dir / "mode_features.tsv.gz", sep="\t", dtype={"subject_id": str})
    feature_columns = [f"{channel}_{mode:03d}" for channel in ["g", "w"] for mode in range(1,130)]
    missing = set(feature_columns) - set(features.columns)
    if missing:
        raise ValueError(f"影像缺少 G/W 特征：{sorted(missing)[:5]}")
    if features.subject_id.duplicated().any():
        raise ValueError("影像特征存在重复 ID")
    summaries = pd.read_csv(feature_dir / "summary.tsv", sep="\t", dtype={"subject_id": str})
    if "lesion_volume_mm3" not in summaries:
        raise ValueError("summary.tsv 必须包含原始掩膜计算的 lesion_volume_mm3，不能用重采样后体积替代")
    if summaries.subject_id.duplicated().any():
        raise ValueError("病灶体积表存在重复 ID")
    if "lesion_volume_mm3" not in features:
        features = features.merge(summaries[["subject_id", "lesion_volume_mm3"]], on="subject_id", how="left", validate="one_to_one")
    numeric = features[feature_columns+["lesion_volume_mm3"]].to_numpy(float)
    if not np.isfinite(numeric).all() or (features.lesion_volume_mm3 <= 0).any():
        raise ValueError("影像特征缺失/非有限或体积非正；失败/空病灶不能以零特征代替")
    if ((features[feature_columns] < -1e-8) | (features[feature_columns] > 1+1e-8)).any().any():
        raise ValueError("G/W 特征超出 [0,1]")
    clinical_ids, imaging_ids = set(mapped.subject_id), set(features.subject_id)
    unmatched = [{"subject_id": value, "reason": "clinical_without_valid_features"} for value in sorted(clinical_ids-imaging_ids)]
    unmatched += [{"subject_id": value, "reason": "features_outside_clinical_cohort"} for value in sorted(imaging_ids-clinical_ids)]
    write_audit(output_dir, "image_clinical_unmatched", unmatched, ["subject_id", "reason"])
    integrated = features.merge(mapped, on="subject_id", how="inner", validate="one_to_one")
    if integrated.empty:
        raise ValueError("临床与有效病灶没有匹配患者")
    integrated = derive_variables(integrated, config["clinical"]["derived"])
    integrated.to_csv(output_dir / "analysis.tsv.gz", sep="\t", index=False, compression="gzip")
    missingness = [{"variable": column, "n": len(integrated), "observed": int(integrated[column].notna().sum()), "missing": int(integrated[column].isna().sum())} for column in integrated if column not in feature_columns]
    write_audit(output_dir, "missingness", missingness)
    definitions = outcome_definitions(config)
    flow = []
    for outcome in definitions:
        if outcome["name"] == "mrs_longitudinal":
            observed = integrated[list(config["statistics"]["families"]["mrs_longitudinal"]["visits"])].notna().any(axis=1)
        else:
            observed = integrated[outcome["name"]].notna()
        flow.append({"outcome": outcome["name"], "family": outcome["family"], "image_matched": len(integrated), "observed_outcome": int(observed.sum()), "missing_outcome": int((~observed).sum()), "covariate_missing": int(integrated.loc[observed,outcome["covariates"]].isna().any(axis=1).sum())})
    write_audit(output_dir, "sample_flow", flow)
    available_variables = [
        column for column in integrated.columns
        if column not in feature_columns and integrated[column].notna().any()
    ]
    result = {"status":"complete", "cohort":cohort, "n":len(integrated), "clinical_n":len(raw),
              "feature_columns":feature_columns, "outcomes":definitions,
              "available_variables": available_variables,
              "variables":config["clinical"]["variables"], "sources":provenance,
              "subject_ids":integrated.subject_id.tolist(), "cohort_subject_ids":raw.subject_id.tolist(),
              "input_csv":"analysis.tsv.gz", "sha256":digest(output_dir / "analysis.tsv.gz")}
    save_json(output_dir / "manifest.json", result)
    return result
