#!/usr/bin/env python3
"""Export aggregate review inputs without copying patient files or changing a run.

Standalone: only NumPy and PyYAML from the project environment are required.
This is an explicit schema exporter, never a recursive directory packer.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml


SCHEMA_VERSION = 1
STAGES = ("clinical", "imputation", "association", "patterns", "validation")
FAMILIES = {"nihss", "mrs_single", "mrs_longitudinal", "moca"}
MONTHS = (3, 6, 12, 24, 36, 48, 60)
NIHSS = tuple("nihss_" + item for item in (
    "1a", "1b", "1c", "2", "3", "4", "5a", "5b", "6a", "6b", "7", "8", "9", "10", "11",
    "consciousness", "motor_left", "motor_right", "motor_total"))
OUTCOMES = set(NIHSS) | {f"mrs_{m}m" for m in MONTHS} | {
    "mrs_longitudinal", "moca_3m", "moca_6m", "moca_12m"}
COVARIATES = {"age", "sex", "log_lesion_volume", "lesion_volume_mm3", "log_volume",
              "prestroke_mrs", "acute_nihss", "education"}
FEATURES = {f"{channel}_{mode:03d}" for channel in ("g", "w") for mode in range(1, 130)}
MODEL_STATUS = {"estimated", "not_estimable", "assumption_limited", "mi_diagnostic_warning",
                "no_image_pattern", "frozen_pending_validation"}
MI_STATUS = {"diagnostic_passed", "no_missing_covariates", "diagnostic_incomplete", "psrf_warning",
             "not_run", "training_mi_warning", "validation_conditional_mi_warning", "not_estimable"}
PO_STATUS = {"not_applicable", "requires_random_effect_threshold_sensitivity", "binary_no_po_constraint",
             "test_unavailable", "assumption_warning", "no_warning_detected"}
ANALYSES = {"gw_primary", "jmax_sensitivity"}
ELIGIBILITY = {"replicated_with_passed_training_and_validation_diagnostics", "training_mi_warning_or_unverified",
               "mi_coefficient_direction_unstable_or_unverified", "validation_model_diagnostic_warning",
               "not_replicated", "direction_not_replicated"}
HEX = re.compile(r"[0-9a-f]{64}\Z")
MISSING = {"", "NA", "NaN", "nan", "null", "None"}
MAX_FILE_BYTES = 128 * 1024 * 1024


class ExportError(ValueError):
    """A safe, value-free diagnostic: never echo rejected clinical content."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ExportError(message)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def contained(path: Path, parent: Path) -> Path:
    resolved = path.resolve()
    require(resolved.is_relative_to(parent.resolve()), "Source path leaves its allowed directory")
    return resolved


def json_load(data: bytes) -> dict:
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "Duplicate JSON key")
            result[key] = value
        return result
    try:
        value = json.loads(data, object_pairs_hook=unique)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ExportError("Invalid JSON metadata") from exc
    require(isinstance(value, dict), "Expected a JSON object")
    return value


def number(value, field: str, *, integer: bool = False, minimum=None, maximum=None, nullable=True):
    if value is None or (isinstance(value, str) and value in MISSING):
        require(nullable, f"Missing required number: {field}")
        return None
    require(not isinstance(value, bool), f"Invalid numerical type: {field}")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ExportError(f"Invalid number: {field}") from exc
    require(math.isfinite(result), f"Non-finite number: {field}")
    require(not integer or result.is_integer(), f"Non-integer count: {field}")
    require(minimum is None or result >= minimum, f"Number below permitted range: {field}")
    require(maximum is None or result <= maximum, f"Number above permitted range: {field}")
    return int(result) if integer else result


def boolean(value, field: str):
    if value is None or (isinstance(value, str) and value in MISSING):
        return None
    if value in (True, 1, "1", "TRUE", "true", "True"):
        return True
    if value in (False, 0, "0", "FALSE", "false", "False"):
        return False
    raise ExportError(f"Invalid boolean: {field}")


def expected_family(name: str) -> str:
    if name in NIHSS:
        return "nihss"
    if name == "mrs_longitudinal":
        return name
    return "mrs_single" if name.startswith("mrs_") else "moca"


def outcome_registry(job: dict) -> dict[str, dict]:
    """Use the frozen job, but accept only this project's canonical semantic names.

    Free-text titles, raw SAS fields, paths, ID rules and raw config are never returned.
    """
    config = job.get("config", {})
    clinical = config.get("clinical", {})
    definitions = job.get("outcomes")
    if definitions is None:
        definitions = []
        for family, spec in config.get("statistics", {}).get("families", {}).items():
            require(family in FAMILIES, "Unknown configured outcome family")
            if not spec.get("enabled", True):
                continue
            covariates = clinical.get("covariate_sets", {}).get(spec.get("covariate_set"), [])
            for name in spec.get("outcomes", []):
                variable = clinical.get("variables", {}).get(name, {})
                definitions.append({"name": name, "family": family,
                    "type": "longitudinal" if name == "mrs_longitudinal" else variable.get("type"),
                    "levels": list(range(7)) if name == "mrs_longitudinal" else variable.get("levels"),
                    "covariates": covariates, "visit_columns": list(spec.get("visits", {})),
                    "visit_times": list(spec.get("visits", {}).values())})
    require(isinstance(definitions, list) and bool(definitions), "Frozen job lacks outcome definitions")
    registry = {}
    for item in definitions:
        require(isinstance(item, dict), "Invalid outcome definition")
        name = item.get("name")
        require(name in OUTCOMES and name not in registry, "Unknown or duplicate canonical outcome")
        family = item.get("family")
        require(family == expected_family(name), "Outcome family is inconsistent")
        kind = item.get("type")
        require(kind in {"continuous", "ordinal", "longitudinal"}, "Unsupported outcome type")
        covariates = item.get("covariates", [])
        require(isinstance(covariates, list) and set(covariates) <= COVARIATES,
                "Unknown canonical covariate")
        levels = item.get("levels")
        if levels is not None:
            require(isinstance(levels, list) and len(levels) <= 64, "Invalid outcome levels")
            levels = [number(x, "levels", integer=True, minimum=0, maximum=100, nullable=False) for x in levels]
        visits = item.get("visit_columns", [])
        require(isinstance(visits, list) and set(visits) <= {f"mrs_{m}m" for m in MONTHS}, "Invalid visit columns")
        times = item.get("visit_times", [])
        require(isinstance(times, list) and len(times) == len(visits), "Invalid visit timing")
        times = [number(x, "visit_time", minimum=0, maximum=120, nullable=False) for x in times]
        registry[name] = {"name": name, "family": family, "type": kind, "levels": levels,
            "covariates": covariates, "visit_columns": visits, "visit_times": times,
            "canonical_direction": "higher_better" if name.startswith("moca_") else "higher_worse",
            "coefficient_unit": "score_points_per_raw_damage_unit" if kind == "continuous" else "log_odds_per_raw_damage_unit",
            "effect_exposure_scale": "unstandardized_g_w_damage"}
    return registry


class ReviewExport:
    def __init__(self, project: Path, config: Path, run: str):
        self.project = project.resolve()
        require(self.project.is_dir(), "Project directory does not exist")
        require(bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", run)) and run not in {".", ".."}, "Unsafe run name")
        self.sources: dict[Path, str] = {}
        self.entries: dict[str, bytes] = {}
        self.files: dict[str, dict] = {}
        self.source_checksums: dict[str, str] = {}
        self.missing: list[dict] = []
        self.sanitized: Counter = Counter()
        self.stage_refs: dict[str, dict] = {}
        self.stages: dict[str, tuple[Path, dict, dict, dict]] = {}
        config_path = contained(config if config.is_absolute() else self.project / config, self.project)
        try:
            settings = yaml.safe_load(self.read(config_path))
        except yaml.YAMLError as exc:
            raise ExportError("Invalid project configuration") from exc
        require(isinstance(settings, dict), "Invalid project configuration object")
        result_base = settings.get("paths", {}).get("results", "results")
        require(isinstance(result_base, str), "Invalid configured results directory")
        self.results = contained(self.project / result_base, self.project)
        require(self.results != self.project and "archive" not in self.results.relative_to(self.project).parts,
                "Results must be a dedicated non-archive project directory")
        self.run_root = contained(self.results / run, self.results)
        require("archive" not in self.run_root.relative_to(self.project).parts, "Archived runs are not export sources")
        self.run_sha = sha(run.encode())
        state = json_load(self.read(contained(self.run_root / "stages.json", self.run_root)))
        require(isinstance(state.get("stages"), dict), "Missing current stage pointers")
        for name in STAGES:
            record = state["stages"].get(name, {})
            require(isinstance(record, dict), "Invalid stage record")
            status = record.get("status", "unavailable")
            status = self.enum(status, {"completed", "failed", "blocked", "pending", "running", "stopped", "stale", "unavailable"}, "stage_status")
            self.stage_refs[name] = {"status": status}
            if status != "completed":
                self.missing.append({"stage": name, "artifact": "stage", "reason": "stage_not_completed"})
                continue
            signature = record.get("signature")
            require(isinstance(signature, str) and HEX.fullmatch(signature), "Invalid completed-stage signature")
            required_path = self.run_root / "stages" / name / signature[:16]
            supplied = record.get("path")
            require(isinstance(supplied, str), "Completed stage lacks an exact path")
            # Existing pipeline paths are absolute; relative paths are project-relative.
            target = contained(self.project / supplied, self.run_root / "stages" / name)
            require(target == required_path.absolute() and target.is_dir(), "Stage pointer is not the current run signature directory")
            complete = json_load(self.read(contained(target / "complete.json", target)))
            require(complete.get("signature") == signature and isinstance(complete.get("outputs"), dict),
                    "Stage completion metadata disagrees with current pointer")
            self.stage_refs[name].update(signature=signature, source_ref=f"stages/{name}/{signature[:16]}")
            self.stages[name] = (target, complete["outputs"], {}, {})
            job_data = self.source(name, "job.json", required=True)
            job = json_load(job_data)
            require(job.get("stage") == name, "Frozen job belongs to another stage")
            registry = outcome_registry(job)
            self.stages[name] = (target, complete["outputs"], job, registry)

    def read(self, path: Path) -> bytes:
        require(path.is_file() and path.stat().st_size <= MAX_FILE_BYTES, "Missing or oversized selected source file")
        data = path.read_bytes()
        digest = sha(data)
        require(path not in self.sources or self.sources[path] == digest, "Source changed while exporting")
        self.sources[path] = digest
        return data

    def enum(self, value, choices: set, field: str) -> str:
        if isinstance(value, str) and value in choices:
            return value
        self.sanitized[field] += 1
        return "unknown"

    def source(self, stage: str, relative: str, *, required=False) -> bytes | None:
        target, outputs, _, _ = self.stages[stage]
        path = contained(target / relative, target)
        if relative not in outputs:
            require(not path.exists(), "Selected artifact is not recorded in completed-stage outputs")
            require(not required, "Required completed-stage artifact is absent")
            self.missing.append({"stage": stage, "artifact": relative, "reason": "artifact_not_produced"})
            return None
        digest = outputs[relative]
        require(isinstance(digest, str) and HEX.fullmatch(digest), "Invalid source artifact checksum")
        data = self.read(path)
        require(sha(data) == digest, "Source artifact checksum mismatch")
        self.source_checksums[f"{stage}/{relative}"] = digest
        return data

    def add(self, name: str, data: bytes, **metadata) -> None:
        require(name not in self.entries and not name.startswith("/") and ".." not in Path(name).parts, "Invalid export entry")
        self.entries[name] = data
        self.files[name] = {"sha256": sha(data), "bytes": len(data), **metadata}

    def add_json(self, name: str, value, **metadata) -> None:
        self.add(name, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False).encode("utf-8") + b"\n", **metadata)

    def add_table(self, name: str, rows: list[dict], columns: list[str]) -> None:
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        self.add(name, buffer.getvalue().encode("utf-8"), rows=len(rows), columns=columns)

    def rows(self, stage: str, relative: str, columns: list[str]) -> list[dict] | None:
        data = self.source(stage, relative)
        if data is None:
            return None
        try:
            reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t" if relative.endswith(".tsv") else ",")
            fields = reader.fieldnames or []
            require(len(fields) == len(set(fields)) and set(columns) <= set(fields), "Selected aggregate table schema is incomplete or duplicated")
            result = []
            for row in reader:
                require(None not in row and all(row.get(key) is not None for key in columns), "Malformed aggregate table row")
                result.append({key: row[key] for key in columns})
            return result
        except (UnicodeError, csv.Error) as exc:
            raise ExportError("Invalid aggregate table encoding") from exc

    def clean(self, row: dict, stage: str) -> dict:
        registry = self.stages[stage][3]
        output = {}
        for key, value in row.items():
            if key in {"outcome", "target", "pattern"}:
                allowed = set(registry)
                if key == "pattern":
                    targets = self.stages[stage][2].get("config", {}).get("patterns", {}).get("targets", [])
                    require(isinstance(targets, list) and set(targets) <= OUTCOMES, "Unknown configured pattern target")
                    allowed = set(targets)
                require(value in allowed, "Aggregate row has an unknown configured outcome")
                output[key] = value
            elif key == "family":
                require(value in FAMILIES and row.get("outcome", row.get("target")) in registry and
                        value == registry[row.get("outcome", row.get("target"))]["family"], "Aggregate row has an inconsistent family")
                output[key] = value
            elif key in {"feature", "term"}:
                require(value in FEATURES or (key == "term" and value == "jmax") or value in MISSING, "Unknown canonical image feature")
                output[key] = None if value in MISSING else value
            elif key == "analysis":
                require(value in ANALYSES, "Unknown association analysis identity")
                output[key] = value
            elif key == "status":
                output[key] = self.enum(value, MODEL_STATUS | MI_STATUS, key)
            elif key == "mi_status":
                output[key] = self.enum(value, MI_STATUS, key)
            elif key == "po_status":
                output[key] = self.enum(value, PO_STATUS, key)
            elif key in {"sparse_category", "selected", "converged_all_folds", "lambda_boundary", "annotation_eligible",
                         "training_diagnostics_passed", "coefficient_direction_stable"}:
                output[key] = boolean(value, key)
            elif key == "annotation_reason":
                parts = value.split(";") if isinstance(value, str) else []
                output[key] = ";".join(parts) if parts and set(parts) <= ELIGIBILITY else self.enum(None, set(), key)
            else:
                integer = key in {"mode", "n", "observations", "planned_family_tests", "imputation", "m", "iterations", "events",
                    "image_matched", "excluded_outcome", "observed_outcome", "missing_outcome", "covariate_missing"}
                positive = integer or key.endswith("standard_error") or key in {"standard_error", "coefficient_sd", "scale", "lambda", "mean_log_loss", "information_condition", "gradient", "maximum_psrf"}
                probability = key in {"q", "p_joint", "p", "po_min_p", "heteroscedasticity_p", "ceiling_fraction", "selection_frequency", "sign_agreement_with_ensemble"} or key.endswith("_p")
                output[key] = number(value, key, integer=integer, minimum=1 if key == "mode" else 0 if positive or probability else None,
                    maximum=129 if key == "mode" else 1 if probability else None, nullable=key != "mode")
        return output

    def table(self, stage: str, relative: str, columns: list[str]) -> None:
        rows = self.rows(stage, relative, columns)
        if rows is not None:
            self.add_table(f"{stage}/{relative.replace('.tsv', '.csv')}", [self.clean(row, stage) for row in rows], columns)

    def aggregates(self) -> None:
        metadata = []
        for name, (_, _, job, registry) in self.stages.items():
            metadata.append({"stage": name, "cohort": "discovery_and_validation" if name == "clinical" else "validation" if name == "validation" else "discovery",
                             "outcomes": list(registry.values()), "scientific_settings": self.scientific_settings(job)})
            if name == "clinical":
                columns = ["outcome", "family", "image_matched", "excluded_outcome", "observed_outcome", "missing_outcome", "covariate_missing"]
                for cohort in ("discovery", "validation"):
                    self.table(name, f"{cohort}/audit/sample_flow.tsv", columns)
            elif name == "imputation":
                records = []
                for outcome in registry:
                    raw = self.source(name, f"{outcome}/diagnostics.json")
                    if raw is None:
                        continue
                    value = json_load(raw)
                    record = {key: value.get(key) for key in ["status", "m", "iterations", "maximum_psrf", "events"]}
                    records.append(self.clean({"outcome": outcome, **record}, name))
                self.add_table("imputation/diagnostics.csv", records, ["outcome", "status", "m", "iterations", "maximum_psrf", "events"])
                self.table(name, "imputation_summary.csv", ["outcome", "status", "n"])
            elif name == "association":
                base = ["outcome", "family", "mode", "analysis", "status", "n", "observations", "p_joint", "mi_status"]
                effects = [f"{channel}_{field}" for channel in ("g", "w") for field in ("estimate", "standard_error", "lower", "upper", "p")]
                self.table(name, "associations.csv", base + effects + ["q", "planned_family_tests"])
                self.table(name, "jmax_sensitivity.csv", base + [f"jmax_{x}" for x in ("estimate", "standard_error", "lower", "upper", "p")] + ["q", "planned_family_tests"])
                self.table(name, "threshold_sensitivity.csv", ["outcome", "family", "mode", "analysis", "threshold", "status", "p_joint", "term", "estimate", "standard_error", "lower", "upper", "p"])
                self.table(name, "time_interaction.csv", ["outcome", "family", "mode", "analysis", "status", "p_joint", "q", "planned_family_tests"])
                self.diagnostics(name)
            elif name == "patterns":
                self.table(name, "patterns.csv", ["target", "status", "lambda", "lambda_boundary"])
                targets = job.get("config", {}).get("patterns", {}).get("targets", [])
                require(isinstance(targets, list) and set(targets) <= set(registry), "Unknown configured pattern target")
                for target in targets:
                    self.table(name, f"{target}/coefficients.csv", ["feature", "mean", "scale", "coefficient", "raw_coefficient", "selection_frequency", "sign_agreement_with_ensemble", "coefficient_sd"])
                    self.table(name, f"{target}/cross_validation.csv", ["lambda", "mean_log_loss", "standard_error", "selected", "converged_all_folds"])
            elif name == "validation":
                columns = ["pattern", "outcome", "family", "status", "n", "p_joint", "estimate", "lower", "upper", "training_diagnostics_passed", "coefficient_direction_stable", "q", "planned_family_tests"]
                self.table(name, "validation.csv", columns)
                self.table(name, "pattern_eligibility.csv", columns + ["annotation_eligible", "annotation_reason"])
                self.table(name, "planned_tests.csv", ["pattern", "outcome", "family"])
        self.add_json("outcomes.json", {"definitions_by_stage": metadata,
            "direction_note": "Canonical score direction; no patient values were inspected or transformed.",
            "comparison_note": "G/W estimates are unstandardized; continuous effects and ordinal log odds are not a common effect scale."})

    def scientific_settings(self, job: dict) -> dict:
        """Select numerical settings and fixed enums, never serialize arbitrary config."""
        config = job.get("config", {})
        groups = {
            "mi": (config.get("statistics", {}).get("mi", {}),
                {"m", "maxit", "initial_iterations", "increment", "iteration_increment", "max_iterations", "maximum_iterations", "seed", "psrf_threshold"},
                {"inference_outcome_predictor", "learning_outcome_predictor"}, {"outcome_policy": {"observed_only"}}),
            "fdr": (config.get("statistics", {}).get("fdr", {}), {"q"}, {"failed_tests_in_denominator"}, {"method": {"BH"}}),
            "ordinal": (config.get("statistics", {}).get("models", {}).get("ordinal", {}),
                {"po_alpha", "sparse_level_minimum", "information_condition_limit", "gradient_tolerance"}, {"threshold_sensitivity"},
                {"package": {"ordinal"}, "link": {"logit"}}),
            "linear": (config.get("statistics", {}).get("models", {}).get("linear", {}), {"ceiling_fraction_warning"}, set(), {"covariance": {"HC3"}}),
            "longitudinal": (config.get("statistics", {}).get("models", {}).get("longitudinal", {}),
                {"n_agq", "optimizer_iterations", "optimizer_relative_tolerance"}, {"enabled", "time_interaction", "threshold_sensitivity"},
                {"time": {"categorical"}, "random_intercept": {"subject_id"}, "optimizer": {"ucminf", "nlminb"}}),
            "patterns": (config.get("patterns", {}), {"folds", "alpha", "lambda_min", "lambda_max", "n_lambda", "seed", "maxiter", "tolerance"},
                {"standardize_within_training", "imputation_uses_outcome", "covariates_penalized", "score_contains_covariates"},
                {"package": {"ordinalNet"}, "selection": {"one_se", "one_standard_error", "minimum"}}),
        }
        result = {}
        for group, (source, numerical, booleans, enums) in groups.items():
            require(isinstance(source, dict), "Invalid scientific settings object")
            selected = {}
            for key in sorted(numerical & set(source)):
                selected[key] = number(source[key], key, minimum=0, nullable=False)
            for key in sorted(booleans & set(source)):
                selected[key] = boolean(source[key], key)
            for key in sorted(set(enums) & set(source)):
                selected[key] = self.enum(source[key], enums[key], "setting_" + key)
            result[group] = selected
        clinical = config.get("clinical", {})
        methods = config.get("statistics", {}).get("mi", {}).get("methods", {})
        require(isinstance(methods, dict) and set(methods) <= {"continuous", "binary", "nominal", "ordinal", "categorical"},
                "Unknown MI variable type")
        result["mi"]["methods"] = {key: self.enum(value, {"pmm", "logreg", "polyreg", "polr"}, "mi_method")
                                   for key, value in sorted(methods.items())}
        derived = clinical.get("derived", {})
        allowed_derived = {"nihss_consciousness", "nihss_motor_left", "nihss_motor_right", "nihss_motor_total", "log_lesion_volume"}
        require(isinstance(derived, dict) and set(derived) <= allowed_derived, "Unknown configured derived variable")
        result["derived"] = {}
        for target, definition in sorted(derived.items()):
            require(isinstance(definition, dict), "Invalid derived variable definition")
            operation = definition.get("operation")
            require(operation == ("log1p" if target == "log_lesion_volume" else "sum_complete"), "Unsupported derived variable operation")
            inputs = definition.get("inputs")
            require(isinstance(inputs, list) and bool(inputs) and len(set(inputs)) == len(inputs) and
                    set(inputs) <= (COVARIATES if target == "log_lesion_volume" else set(NIHSS)),
                    "Unknown canonical derived input")
            result["derived"][target] = {"operation": operation, "inputs": inputs}
        result["death_rules"] = {}
        for cohort in ("discovery", "validation"):
            effective = {**clinical.get("death", {}), **clinical.get("cohorts", {}).get(cohort, {}).get("death", {})}
            choices = {"policy": {"mrs_equals_6"}, "moca_conflict_policy": {"error", "structural_missing"},
                       "mrs_conflict_policy": {"error", "exclude_mrs_moca"}}
            result["death_rules"][cohort] = {key: self.enum(effective[key], values, "setting_" + key)
                for key, values in choices.items() if key in effective}
        return result

    def diagnostics(self, stage: str) -> None:
        keys = ["outcome", "family", "mode", "analysis"]
        columns = keys + ["imputation", "po_status", "po_min_p", "sparse_category", "information_condition", "gradient", "heteroscedasticity_p", "ceiling_fraction"]
        rows = self.rows(stage, "model_diagnostics.csv", columns)
        if rows is None:
            return
        groups = {}
        for raw in rows:
            row = self.clean(raw, stage)
            identity = tuple(row[x] for x in keys)
            item = groups.setdefault(identity, {**{x: row[x] for x in keys}, "fits": 0,
                **{"po_" + value: 0 for value in sorted(PO_STATUS | {"unknown"})}, "sparse_fits": 0,
                "max_information_condition": None, "max_gradient": None, "ceiling_fraction_max": None,
                "heteroscedasticity_p_min": None, "po_min_p": None})
            item["fits"] += 1
            item["po_" + row["po_status"]] += 1
            item["sparse_fits"] += int(row["sparse_category"] is True)
            for source, target, operation in [("information_condition", "max_information_condition", max),
                ("gradient", "max_gradient", max), ("ceiling_fraction", "ceiling_fraction_max", max),
                ("heteroscedasticity_p", "heteroscedasticity_p_min", min), ("po_min_p", "po_min_p", min)]:
                value = row[source]
                if value is not None:
                    item[target] = value if item[target] is None else operation(value, item[target])
        output = list(groups.values())
        cols = keys + ["fits"] + ["po_" + value for value in sorted(PO_STATUS | {"unknown"})] + ["sparse_fits", "max_information_condition", "max_gradient", "ceiling_fraction_max", "heteroscedasticity_p_min", "po_min_p"]
        self.add_table("association/model_diagnostics_summary.csv", output, cols)

    def reference(self) -> str:
        directory = contained(self.project / "reference", self.project)
        manifest = json_load(self.read(contained(directory / "manifest.json", directory)))
        reference_id = manifest.get("reference_id")
        require(isinstance(reference_id, str) and HEX.fullmatch(reference_id), "Invalid healthy reference identity")
        files = manifest.get("files", {})

        def reference_file(name):
            value = self.read(contained(directory / name, directory))
            require(files.get(name) == sha(value), "Healthy reference file checksum mismatch")
            return value

        for stage, (_, _, job, _) in self.stages.items():
            if stage != "clinical":
                require(job.get("reference_id") == reference_id, "Stage and healthy reference identity differ")
        try:
            with np.load(io.BytesIO(reference_file("modes.npz")), allow_pickle=False) as archive:
                values, vectors, active, positive, strength = [archive[key] for key in ("values", "vectors", "active", "positive", "strength")]
            with np.load(io.BytesIO(reference_file("sc.npz")), allow_pickle=False) as archive:
                weights = archive["weights"]
        except (KeyError, ValueError, OSError) as exc:
            raise ExportError("Invalid healthy numeric reference archive") from exc
        require(values.shape == (130,) and vectors.shape == (135, 130) and active.shape == (135,) and positive.shape == (130,) and strength.shape == (135,) and weights.shape == (135, 135), "Healthy reference dimensions differ from LQT135")
        require(active.dtype.kind == "b" and positive.dtype.kind == "b" and int(active.sum()) == 130 and int(positive.sum()) == 129, "Invalid healthy reference masks")
        require(all(array.dtype.kind in "fi" and np.isfinite(array).all() for array in (values, vectors, strength, weights)), "Non-numerical healthy reference")
        require(np.array_equal(np.flatnonzero(positive), np.arange(1, 130)) and np.all(np.diff(values) >= -1e-12) and np.all(values[positive] > 0), "Invalid frozen spectral ordering")
        require(np.all(weights >= 0) and np.allclose(weights, weights.T, atol=1e-12) and np.all(np.diag(weights) == 0), "Invalid healthy SC weights")
        require(np.allclose(strength, weights.sum(axis=1)) and np.array_equal(active, strength > 0) and np.count_nonzero(np.triu(weights > 0, 1)) == 1375, "SC and spectral support disagree")
        selected = vectors[:, positive]
        require(np.allclose(selected.T @ selected, np.eye(129), atol=1e-6), "Healthy modes are not orthonormal")
        scale = np.divide(1., np.sqrt(strength), out=np.zeros_like(strength), where=strength > 0)
        laplacian = np.diag(active.astype(float)) - scale[:, None] * weights * scale[None, :]
        require(np.max(np.abs(laplacian @ selected - selected * values[positive])) < 1e-6, "Healthy modes do not match SC weights")
        buffer = io.BytesIO()
        np.savez_compressed(buffer, mode_ids=np.arange(1, 130), values=values[positive], vectors=selected,
                            active=active, strength=strength, weights=weights)
        self.add("reference/healthy_modes.npz", buffer.getvalue(), arrays={"mode_ids": [129], "values": [129], "vectors": [135, 129], "active": [135], "strength": [135], "weights": [135, 135]})
        node_rows = list(csv.DictReader(io.StringIO(reference_file("nodes.tsv").decode("utf-8")), delimiter="\t"))
        labels = canonical_nodes()
        require(len(node_rows) == len(labels) == 135, "Invalid healthy node count")
        clean_nodes = []
        for index, (row, label) in enumerate(zip(node_rows, labels), 1):
            require(row.get("label") == label and number(row.get("index"), "node_index", integer=True, nullable=False) == index, "Healthy node labels/order differ from the supported atlas")
            cleaned = {"index": index, "label": label, "active": bool(active[index - 1])}
            for key, choices in {"group": {"cortex", "subcortex", "cerebellum", "brainstem"},
                "atlas_name": {"Schaefer100", "AAL", "HarvardOxford"},
                "coverage_type": {"gray_parcel_coverage", "mixed_tissue_region_coverage"}}.items():
                require(row.get(key) in choices, "Unknown healthy node category")
                cleaned[key] = row[key]
            cleaned["network_id"] = number(row.get("network_id"), "network_id", integer=True, minimum=1, maximum=11, nullable=False)
            for key in ("x_mm", "y_mm", "z_mm"):
                cleaned[key] = number(row.get(key), key, minimum=-200, maximum=200, nullable=False)
            clean_nodes.append(cleaned)
        self.add_table("reference/nodes.csv", clean_nodes, list(clean_nodes[0]))
        # Rebuild summaries from the validated healthy arrays, not arbitrary strings/extra NPZ members.
        spectrum = []
        for position in range(129):
            squared = selected[:, position] ** 2
            record = {"mode": position + 1, "eigenvalue": float(values[positive][position]), "effective_nodes": float(1 / np.sum(squared ** 2))}
            for group in ("cortex", "subcortex", "cerebellum", "brainstem"):
                record[group + "_loading"] = float(sum(squared[i] for i, row in enumerate(clean_nodes) if row["group"] == group))
            spectrum.append(record)
        self.add_table("reference/spectrum.csv", spectrum, list(spectrum[0]))
        return reference_id

    def build(self) -> dict:
        require(bool(self.stages), "No completed clinical/statistical stage is available")
        self.aggregates()
        reference_id = self.reference()
        manifest = {"schema_version": SCHEMA_VERSION, "exported_at_utc": datetime.now(timezone.utc).isoformat(),
            "run_name_sha256": self.run_sha, "reference_id": reference_id, "stages": self.stage_refs,
            "files": self.files, "verified_source_sha256": self.source_checksums,
            "missing_artifacts": self.missing, "sanitized_enum_counts": dict(self.sanitized),
            "privacy_contract": "aggregate_only_no_patient_rows_no_source_config_no_raw_rds",
            "diagnostic_note": "Minimum diagnostic p values describe warnings across fits; they are not pooled hypothesis tests.",
            "verification_scope": "Current completed-stage pointers and selected source checksums; no freshness inference or model refit."}
        # Manifest cannot contain its own checksum; every data entry is listed above.
        self.entries["manifest.json"] = json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False).encode() + b"\n"
        return manifest

    def verify_unchanged(self) -> None:
        for path, digest in self.sources.items():
            require(path.is_file() and sha(path.read_bytes()) == digest, "Source changed during export; retry after a stable snapshot")

    def validate_output(self, output: Path) -> Path:
        output = output.absolute()
        require(output.suffix.lower() == ".zip" and output.parent.is_dir(), "Output must be a ZIP in an existing directory")
        require(not output.is_symlink() and not output.exists(), "Output already exists; choose a new ZIP name")
        resolved = output.resolve()
        require(not resolved.is_relative_to(self.run_root) and not resolved.is_relative_to(self.project / "reference"), "Output cannot be inside a source run or reference directory")
        return output

    def save(self, output: Path) -> None:
        output = self.validate_output(output)
        self.verify_unchanged()
        # Exclusive creation prevents overwrite even if another process wins the race.
        created_identity = None
        try:
            with output.open("xb") as handle:
                stat = os.fstat(handle.fileno())
                created_identity = (stat.st_dev, stat.st_ino)
                with zipfile.ZipFile(handle, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
                    for name, data in self.entries.items():
                        archive.writestr(name, data)
            self.verify_unchanged()
        except FileExistsError as exc:
            raise ExportError("Output already exists; choose a new ZIP name") from exc
        except BaseException:
            if created_identity is not None and output.exists():
                current = output.lstat()
                if (current.st_dev, current.st_ino) == created_identity:
                    output.unlink()
            raise


def canonical_nodes() -> list[str]:
    """Fixed public atlas labels; a patient identifier cannot become an exported node name."""
    hemispheres = {
        "LH": [("Vis", 9), ("SomMot", 6), ("DorsAttn_Post", 6), ("DorsAttn_PrCv", 1), ("DorsAttn_FEF", 1),
            ("SalVentAttn_ParOper", 1), ("SalVentAttn_FrOperIns", 2), ("SalVentAttn_PFCl", 1), ("SalVentAttn_Med", 3),
            ("Limbic_OFC", 1), ("Limbic_TempPole", 2), ("Cont_Par", 1), ("Cont_PFCl", 1), ("Cont_pCun", 1), ("Cont_Cing", 1),
            ("Default_Temp", 2), ("Default_Par", 2), ("Default_PFC", 7), ("Default_pCunPCC", 2)],
        "RH": [("Vis", 8), ("SomMot", 8), ("DorsAttn_Post", 5), ("DorsAttn_PrCv", 1), ("DorsAttn_FEF", 1),
            ("SalVentAttn_TempOccPar", 2), ("SalVentAttn_FrOperIns", 1), ("SalVentAttn_Med", 2),
            ("Limbic_OFC", 1), ("Limbic_TempPole", 1), ("Cont_Par", 2), ("Cont_PFCl", 4), ("Cont_Cing", 1), ("Cont_PFCmp", 1), ("Cont_pCun", 1),
            ("Default_Par", 1), ("Default_Temp", 3), ("Default_PFCv", 2), ("Default_PFCdPFCm", 3), ("Default_pCunPCC", 2)]}
    labels = [f"7Networks_{side}_{region}_{index}" for side, specs in hemispheres.items()
              for region, count in specs for index in range(1, count + 1)]
    labels += [f"{side}_{region}" for region in ("Thalamus", "Caudate", "Putamen", "Pallidum") for side in ("LH", "RH")]
    labels += [f"Cerebelum_{region}_{side}" for region in ("Crus1", "Crus2", "3", "4_5", "6", "7b", "8", "9", "10") for side in ("L", "R")]
    labels += [f"Vermis_{region}" for region in ("1_2", "3", "4_5", "6", "7", "8", "9", "10")] + ["Brainstem"]
    return labels


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config.workstation.yaml"))
    parser.add_argument("--run", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true", help="Validate and preview the bundle without writing anything")
    args = parser.parse_args(argv)
    try:
        exporter = ReviewExport(args.project, args.config, args.run)
        manifest = exporter.build()
        exporter.verify_unchanged()
        exporter.validate_output(args.output)
        if not args.dry_run:
            exporter.save(args.output)
        summary = {"status": "preview" if args.dry_run else "exported", "file_count": len(exporter.entries),
            "reference_id": manifest["reference_id"], "sanitized_enum_counts": manifest["sanitized_enum_counts"]}
        if args.dry_run:
            summary.update(uncompressed_bytes=sum(map(len, exporter.entries.values())), files=manifest["files"],
                           missing_artifacts=manifest["missing_artifacts"])
        else:
            summary.update(output=str(args.output.resolve()), zip_bytes=args.output.stat().st_size,
                missing_artifact_count=len(manifest["missing_artifacts"]),
                missing_by_reason=dict(Counter(item["reason"] for item in manifest["missing_artifacts"])))
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0
    except (ExportError, OSError, TypeError, KeyError) as exc:
        # OSError messages may contain private source paths; only our explicit safe text is shown.
        message = str(exc) if isinstance(exc, ExportError) else type(exc).__name__
        print("Export refused: " + message, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
