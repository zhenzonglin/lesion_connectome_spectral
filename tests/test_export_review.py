"""Synthetic-only privacy, provenance and integrity tests for the review exporter."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pytest

from export_review import ExportError, ReviewExport, canonical_nodes, main


SECRET = "SECRET_PATIENT_C001P00016_PRIVATE_PATH"
RUN = "main_synthetic"
REFERENCE_ID = "f" * 64


def encoded(value):
    return json.dumps(value).encode()


def table_bytes(rows, columns=None, delimiter=","):
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns or list(rows[0]), delimiter=delimiter)
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


@pytest.fixture(scope="module")
def reference_files():
    # A deterministic connected 130-node synthetic graph with the production shapes.
    weights = np.zeros((135, 135))
    edges = {(min(i, (i + 1) % 130), max(i, (i + 1) % 130)) for i in range(130)}
    for i in range(130):
        for j in range(i + 1, 130):
            if len(edges) < 1375:
                edges.add((i, j))
    for i, j in edges:
        weights[i, j] = weights[j, i] = 1 + (i + j) / 1000
    strength = weights.sum(axis=1)
    active = strength > 0
    scale = np.divide(1., np.sqrt(strength), out=np.zeros_like(strength), where=active)
    laplacian = np.diag(active.astype(float)) - scale[:, None] * weights * scale[None, :]
    values, basis = np.linalg.eigh(laplacian[:130, :130])
    vectors = np.zeros((135, 130))
    vectors[:130] = basis
    positive = np.arange(130) > 0
    modes = io.BytesIO()
    np.savez_compressed(modes, values=values, vectors=vectors, active=active, positive=positive,
        strength=strength, patient_ids=np.array([SECRET], dtype=object))
    graph = io.BytesIO()
    np.savez_compressed(graph, weights=weights, patient_ids=np.array([SECRET], dtype=object))
    rows = []
    for index, label in enumerate(canonical_nodes(), 1):
        group = "cortex" if index <= 100 else "subcortex" if index <= 108 else "cerebellum" if index < 135 else "brainstem"
        rows.append({"index": index, "label": label, "group": group,
            "atlas_name": "Schaefer100" if index <= 100 else "AAL" if index < 135 else "HarvardOxford",
            "coverage_type": "mixed_tissue_region_coverage" if index == 135 else "gray_parcel_coverage",
            "network_id": 1 if index <= 100 else 8 if index <= 108 else 10 if index < 135 else 11,
            "x_mm": index % 30, "y_mm": index % 40, "z_mm": index % 50, "private_extra": SECRET})
    return {"modes.npz": modes.getvalue(), "sc.npz": graph.getvalue(), "nodes.tsv": table_bytes(rows, delimiter="\t")}


def outcomes():
    return [
        {"name": "mrs_3m", "family": "mrs_single", "type": "ordinal", "levels": list(range(7)), "covariates": ["age", "sex"], "private": SECRET},
        {"name": "moca_3m", "family": "moca", "type": "continuous", "levels": None, "covariates": ["age", "education"]},
    ]


def assoc_row(mode=1, status="estimated"):
    result = {"outcome": "mrs_3m", "family": "mrs_single", "mode": mode, "analysis": "gw_primary",
        "status": status, "n": 100, "observations": 100, "p_joint": .01 if status == "estimated" else "",
        "mi_status": "diagnostic_passed", "q": .03 if status == "estimated" else "", "planned_family_tests": 129,
        "subject_id": SECRET, "reason": SECRET, "warnings": SECRET}
    for channel in ("g", "w"):
        for field in ("estimate", "standard_error", "lower", "upper", "p"):
            result[channel + "_" + field] = .05 if status == "estimated" else ""
    return result


def add_stage(project, name, files, *, status="completed", job=None):
    signature = hashlib.sha256(name.encode()).hexdigest()
    folder = project / "results" / RUN / "stages" / name / signature[:16]
    folder.mkdir(parents=True)
    if job is None:
        job = {"stage": name, "outcomes": outcomes(), "reference_id": REFERENCE_ID,
            "config": {"patterns": {"targets": ["mrs_3m"]}, "clinical": {"variables": {"mrs_3m": {"title": SECRET, "field": SECRET}}}},
            "input_csv": "/private/" + SECRET, "output_dir": str(folder), "subject_ids": [SECRET]}
    files = {"job.json": encoded(job), **files}
    for relative, data in files.items():
        path = folder / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    outputs = {relative: hashlib.sha256(data).hexdigest() for relative, data in files.items()}
    (folder / "complete.json").write_bytes(encoded({"signature": signature, "outputs": outputs}))
    state_path = project / "results" / RUN / "stages.json"
    state = json.loads(state_path.read_bytes()) if state_path.exists() else {"stages": {}}
    state["stages"][name] = {"signature": signature, "path": str(folder), "status": status, "private": SECRET}
    state_path.write_bytes(encoded(state))
    return folder


@pytest.fixture
def project(tmp_path, reference_files):
    root = tmp_path / "project"
    root.mkdir()
    (root / "config.workstation.yaml").write_text("paths:\n  results: results\nsecret: " + SECRET)
    ref = root / "reference"
    ref.mkdir()
    for name, data in reference_files.items():
        (ref / name).write_bytes(data)
    (ref / "manifest.json").write_bytes(encoded({"reference_id": REFERENCE_ID,
        "files": {name: hashlib.sha256(data).hexdigest() for name, data in reference_files.items()}, "private": SECRET}))
    diagnostics = []
    for index, po in [(1, "no_warning_detected"), (2, "assumption_warning"), (3, SECRET)]:
        diagnostics.append({"outcome": "mrs_3m", "family": "mrs_single", "mode": 1, "analysis": "gw_primary",
            "imputation": index, "po_status": po, "po_min_p": .01 * index, "sparse_category": index == 2,
            "information_condition": 100 * index, "gradient": .000001 * index,
            "heteroscedasticity_p": "", "ceiling_fraction": "", "warnings": SECRET, "subject_id": SECRET})
    add_stage(root, "association", {"associations.csv": table_bytes([assoc_row(), assoc_row(2, "not_estimable")]),
        "model_diagnostics.csv": table_bytes(diagnostics), "private.rds": SECRET.encode(), "report.html": SECRET.encode()})
    return root


def build(project):
    exporter = ReviewExport(project, Path("config.workstation.yaml"), RUN)
    manifest = exporter.build()
    return exporter, manifest


def rewrite_source(project, stage, name, data):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"][stage]["path"])
    (folder / name).write_bytes(data)
    complete = json.loads((folder / "complete.json").read_text())
    complete["outputs"][name] = hashlib.sha256(data).hexdigest()
    (folder / "complete.json").write_bytes(encoded(complete))


def read_csv(data):
    return list(csv.DictReader(io.StringIO(data.decode())))


def test_safe_whitelist_keeps_failed_rows_and_rebuilds_npz(project):
    before = {p: p.read_bytes() for p in project.rglob("*") if p.is_file()}
    exporter, manifest = build(project)
    destination = project / "review.zip"
    exporter.save(destination)
    assert all(p.read_bytes() == data for p, data in before.items())
    with zipfile.ZipFile(destination) as archive:
        assert all(SECRET.encode() not in archive.read(name) for name in archive.namelist())
        assert not any(name.endswith((".rds", ".html")) for name in archive.namelist())
        rows = read_csv(archive.read("association/associations.csv"))
        assert len(rows) == 2 and rows[1]["status"] == "not_estimable" and rows[1]["p_joint"] == ""
        assert "subject_id" not in rows[0] and "reason" not in rows[0]
        with np.load(io.BytesIO(archive.read("reference/healthy_modes.npz")), allow_pickle=False) as arrays:
            assert set(arrays.files) == {"mode_ids", "values", "vectors", "active", "strength", "weights"}
            assert arrays["vectors"].shape == (135, 129)
            assert np.array_equal(arrays["mode_ids"], np.arange(1, 130))
    assert manifest["files"]["association/associations.csv"]["rows"] == 2
    assert manifest["stages"]["association"]["source_ref"].startswith("stages/association/")
    assert str(project) not in json.dumps(manifest)


def test_diagnostics_are_aggregated_and_unknown_enums_redacted(project):
    exporter, manifest = build(project)
    row = read_csv(exporter.entries["association/model_diagnostics_summary.csv"])[0]
    assert row["fits"] == "3" and row["sparse_fits"] == "1" and row["po_unknown"] == "1"
    assert float(row["max_information_condition"]) == 300
    assert manifest["sanitized_enum_counts"]["po_status"] == 1
    assert SECRET.encode() not in exporter.entries["outcomes.json"]


def test_unknown_model_status_redacted_without_dropping_row(project):
    rewrite_source(project, "association", "associations.csv", table_bytes([assoc_row(status=SECRET)]))
    exporter, _ = build(project)
    assert read_csv(exporter.entries["association/associations.csv"])[0]["status"] == "unknown"


@pytest.mark.parametrize("field,value", [("outcome", SECRET), ("family", SECRET), ("analysis", SECRET),
    ("mode", 0), ("mode", 130), ("g_estimate", SECRET), ("g_standard_error", -1), ("p_joint", 1.01)])
def test_invalid_allowed_cells_are_not_exported(project, field, value):
    row = assoc_row()
    row[field] = value
    rewrite_source(project, "association", "associations.csv", table_bytes([row]))
    with pytest.raises(ExportError) as error:
        build(project)
    assert SECRET not in str(error.value)


def test_missing_required_column_is_refused(project):
    row = assoc_row()
    del row["g_standard_error"]
    rewrite_source(project, "association", "associations.csv", table_bytes([row]))
    with pytest.raises(ExportError, match="schema"):
        build(project)


def test_source_hash_mismatch_is_refused(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    (Path(state["stages"]["association"]["path"]) / "associations.csv").write_text(SECRET)
    with pytest.raises(ExportError, match="checksum"):
        build(project)


@pytest.mark.parametrize("variant", ["other_run", "archive", "symlink"])
def test_stage_path_escape_is_refused(project, tmp_path, variant):
    path = project / "results" / RUN / "stages.json"
    state = json.loads(path.read_text())
    old = Path(state["stages"]["association"]["path"])
    if variant == "symlink":
        target = tmp_path / "outside"
        old.rename(target)
        old.symlink_to(target, target_is_directory=True)
    else:
        state["stages"]["association"]["path"] = str(project / ("archive" if variant == "archive" else "results/other_run") / "stage")
        path.write_bytes(encoded(state))
    with pytest.raises(ExportError):
        build(project)


def test_optional_artifact_symlink_escape_is_refused(project, tmp_path):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"]["association"]["path"])
    external = tmp_path / "private.csv"
    external.write_text(SECRET)
    (folder / "jmax_sensitivity.csv").symlink_to(external)
    with pytest.raises(ExportError, match="leaves"):
        build(project)


def test_no_overwrite_and_no_zip_inside_source_run(project):
    exporter, _ = build(project)
    output = project / "exists.zip"
    output.write_bytes(b"preserve")
    with pytest.raises(ExportError, match="already exists"):
        exporter.save(output)
    assert output.read_bytes() == b"preserve"
    with pytest.raises(ExportError, match="source run"):
        exporter.save(project / "results" / RUN / "new.zip")


def test_changed_source_refused_and_partial_zip_removed(project, monkeypatch):
    exporter, _ = build(project)
    original = exporter.verify_unchanged
    calls = []
    def verify():
        calls.append(1)
        if len(calls) == 2:
            raise ExportError("Source changed during export")
        original()
    monkeypatch.setattr(exporter, "verify_unchanged", verify)
    output = project / "changed.zip"
    with pytest.raises(ExportError, match="changed"):
        exporter.save(output)
    assert not output.exists()


def test_dry_run_is_read_only(project, capsys):
    before = {p: p.read_bytes() for p in project.rglob("*") if p.is_file()}
    output = project / "preview.zip"
    assert main(["--project", str(project), "--run", RUN, "--output", str(output), "--dry-run"]) == 0
    assert not output.exists()
    assert all(p.read_bytes() == value for p, value in before.items())
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "preview" and SECRET not in json.dumps(result)


@pytest.mark.parametrize("variant", ["wrong_reference", "tampered_weights", "private_node_label"])
def test_reference_mismatch_or_malformed_reference_refused(project, variant):
    ref = project / "reference"
    manifest = json.loads((ref / "manifest.json").read_text())
    if variant == "wrong_reference":
        manifest["reference_id"] = "e" * 64
    elif variant == "tampered_weights":
        np.savez_compressed(ref / "sc.npz", weights=np.zeros((135, 135)))
        manifest["files"]["sc.npz"] = hashlib.sha256((ref / "sc.npz").read_bytes()).hexdigest()
    else:
        path = ref / "nodes.tsv"
        path.write_text(path.read_text().replace(canonical_nodes()[0], SECRET))
        manifest["files"]["nodes.tsv"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (ref / "manifest.json").write_bytes(encoded(manifest))
    with pytest.raises(ExportError) as error:
        build(project)
    assert SECRET not in str(error.value)


def test_noncompleted_stage_is_omitted_without_archive_fallback(project):
    add_stage(project, "imputation", {"imputation_summary.csv": table_bytes([{"outcome": "mrs_3m", "status": "diagnostic_passed", "n": 100}])}, status="running")
    exporter, manifest = build(project)
    assert not any(name.startswith("imputation/") for name in exporter.entries)
    assert {"stage": "imputation", "artifact": "stage", "reason": "stage_not_completed"} in manifest["missing_artifacts"]


def test_aggregate_clinical_mi_pattern_validation_schemas(project):
    flow = {"outcome": "mrs_3m", "family": "mrs_single", "image_matched": 100, "excluded_outcome": 2,
        "observed_outcome": 80, "missing_outcome": 18, "covariate_missing": 5, "subject_id": SECRET}
    add_stage(project, "clinical", {"discovery/audit/sample_flow.tsv": table_bytes([flow], delimiter="\t"),
        "discovery/audit/death_conflicts.tsv": SECRET.encode(), "discovery/manifest.json": encoded({"subject_ids": [SECRET]})})
    add_stage(project, "imputation", {"mrs_3m/diagnostics.json": encoded({"status": "diagnostic_passed", "m": 30, "iterations": 20,
        "maximum_psrf": 1.01, "events": 0, "subject_ids": [SECRET]}), "mrs_3m/imputations.rds": SECRET.encode(),
        "imputation_summary.csv": table_bytes([{"outcome": "mrs_3m", "status": "diagnostic_passed", "n": 80, "reason": SECRET}])})
    coefficients = {"feature": "g_001", "mean": .01, "scale": .2, "coefficient": .4, "raw_coefficient": 2,
        "selection_frequency": .8, "sign_agreement_with_ensemble": .9, "coefficient_sd": .2, "private": SECRET}
    cv = {"lambda": .1, "mean_log_loss": .6, "standard_error": .1, "selected": True, "converged_all_folds": True}
    add_stage(project, "patterns", {"patterns.csv": table_bytes([{"target": "mrs_3m", "status": "frozen_pending_validation", "lambda": .1, "lambda_boundary": False, "reason": SECRET}]),
        "mrs_3m/coefficients.csv": table_bytes([coefficients]), "mrs_3m/cross_validation.csv": table_bytes([cv]),
        "mrs_3m/folds.csv": SECRET.encode(), "mrs_3m/discovery_scores.csv": SECRET.encode()})
    valid = {"pattern": "mrs_3m", "outcome": "moca_3m", "family": "moca", "status": "estimated", "n": 80,
        "p_joint": .02, "estimate": -.3, "lower": -.5, "upper": -.1, "training_diagnostics_passed": True,
        "coefficient_direction_stable": True, "q": .04, "planned_family_tests": 3, "reason": SECRET}
    add_stage(project, "validation", {"validation.csv": table_bytes([valid]),
        "planned_tests.csv": table_bytes([{key: valid[key] for key in ("pattern", "outcome", "family")}]),
        "pattern_eligibility.csv": table_bytes([{**valid, "annotation_eligible": False, "annotation_reason": "direction_not_replicated"}])})
    exporter, _ = build(project)
    assert "clinical/discovery/audit/sample_flow.csv" in exporter.entries
    assert "imputation/diagnostics.csv" in exporter.entries
    assert "patterns/mrs_3m/coefficients.csv" in exporter.entries
    assert "validation/validation.csv" in exporter.entries
    assert all(SECRET.encode() not in data for data in exporter.entries.values())


def test_unknown_configured_outcome_is_refused(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"]["association"]["path"])
    job = json.loads((folder / "job.json").read_text())
    job["outcomes"][0]["name"] = SECRET
    rewrite_source(project, "association", "job.json", encoded(job))
    with pytest.raises(ExportError, match="canonical outcome"):
        build(project)


def test_scientific_metadata_uses_frozen_job_not_current_config(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"]["association"]["path"])
    job = json.loads((folder / "job.json").read_text())
    job["config"]["statistics"] = {"mi": {"m": 30, "seed": 123, "outcome_policy": SECRET, "private": SECRET},
        "fdr": {"method": "BH", "q": .05, "failed_tests_in_denominator": True}}
    rewrite_source(project, "association", "job.json", encoded(job))
    with (project / "config.workstation.yaml").open("a") as handle:
        handle.write("\nstatistics:\n  mi:\n    m: 999\n")
    exporter, manifest = build(project)
    metadata = json.loads(exporter.entries["outcomes.json"])["definitions_by_stage"][0]
    assert metadata["scientific_settings"]["mi"]["m"] == 30
    assert metadata["scientific_settings"]["mi"]["outcome_policy"] == "unknown"
    assert manifest["sanitized_enum_counts"]["setting_outcome_policy"] == 1
    assert SECRET.encode() not in exporter.entries["outcomes.json"]


def test_completion_signature_mismatch_is_refused(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    path = Path(state["stages"]["association"]["path"]) / "complete.json"
    complete = json.loads(path.read_text())
    complete["signature"] = "0" * 64
    path.write_bytes(encoded(complete))
    with pytest.raises(ExportError, match="completion metadata"):
        build(project)


def test_untracked_selected_artifact_is_refused(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"]["association"]["path"])
    (folder / "jmax_sensitivity.csv").write_text(SECRET)
    with pytest.raises(ExportError, match="not recorded"):
        build(project)


def test_output_symlink_is_never_followed(project, tmp_path):
    exporter, _ = build(project)
    private = tmp_path / "keep.zip"
    private.write_bytes(b"keep")
    output = project / "review.zip"
    output.symlink_to(private)
    with pytest.raises(ExportError, match="already exists"):
        exporter.save(output)
    assert private.read_bytes() == b"keep"


def test_standalone_cli_runs_outside_repository(project, tmp_path):
    output = project / "outside_preview.zip"
    result = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve().parents[1] / "export_review.py"),
        "--project", str(project), "--config", "config.workstation.yaml", "--run", RUN,
        "--output", str(output), "--dry-run"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["status"] == "preview"
    assert SECRET not in result.stdout + result.stderr and not output.exists()


def test_frozen_derived_formulas_and_mi_methods_are_allowlisted(project):
    state = json.loads((project / "results" / RUN / "stages.json").read_text())
    folder = Path(state["stages"]["association"]["path"])
    job = json.loads((folder / "job.json").read_text())
    job["config"]["clinical"]["derived"] = {"nihss_motor_left": {"operation": "sum_complete", "inputs": ["nihss_5a", "nihss_6a"], "private": SECRET}}
    job["config"]["statistics"] = {"mi": {"methods": {"ordinal": "polr", "continuous": SECRET}}}
    rewrite_source(project, "association", "job.json", encoded(job))
    exporter, manifest = build(project)
    metadata = json.loads(exporter.entries["outcomes.json"])["definitions_by_stage"][0]["scientific_settings"]
    assert metadata["derived"]["nihss_motor_left"] == {"operation": "sum_complete", "inputs": ["nihss_5a", "nihss_6a"]}
    assert metadata["mi"]["methods"] == {"ordinal": "polr", "continuous": "unknown"}
    assert manifest["sanitized_enum_counts"]["mi_method"] == 1
    assert SECRET.encode() not in exporter.entries["outcomes.json"]
    job["config"]["clinical"]["derived"]["nihss_motor_left"]["inputs"] = [SECRET]
    rewrite_source(project, "association", "job.json", encoded(job))
    with pytest.raises(ExportError, match="derived input"):
        build(project)


def test_successful_cli_prints_compact_destination_summary(project, capsys):
    output = project / "review.zip"
    assert main(["--project", str(project), "--run", RUN, "--output", str(output)]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "exported" and summary["output"] == str(output.resolve())
    assert summary["zip_bytes"] == output.stat().st_size and "files" not in summary
    assert summary["missing_artifact_count"] > 0
