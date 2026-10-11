"""临床数据工程测试；全部为合成记录，不形成临床结果。"""

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import pandas as pd
import pytest
from functions.configuration import load_config
from functions import clinical


@pytest.fixture
def config():
    root = Path(__file__).resolve().parents[1]
    value = load_config(root / "config.yaml")
    for family in value["statistics"]["families"].values():
        family["enabled"] = False
    value["statistics"]["families"]["nihss"]["enabled"] = True
    value["patterns"]["enabled"] = False
    value["clinical"]["cohorts"]["discovery"].update({"confirmed": True, "rules_confirmed":{"id_mapping":True,"mrs_death":True}, "id":{"field":"id","kind":"string"}})
    return value


def raw_frame(config):
    frame = pd.DataFrame({"id":["0001","a002","0003"]})
    for name, specification in config["clinical"]["variables"].items():
        if "field" in specification:
            frame[specification["field"]] = 0 if "levels" not in specification else specification["levels"][0]
    frame["AGE"] = [50, np.nan, 70]
    return frame


def make_root(tmp_path):
    (tmp_path / "clinical/discovery").mkdir(parents=True)
    (tmp_path / "clinical/discovery/data.sas7bdat").write_bytes(b"synthetic_reader_fixture")


def test_default_contract():
    value = load_config(Path(__file__).resolve().parents[1] / "config.yaml")
    assert len(value["patterns"]["targets"]) == 20
    assert len(value["statistics"]["families"]) == 4
    assert value["clinical"]["covariate_sets"]["moca"] == ["age","sex","log_lesion_volume","education"]
    assert "prestroke_mrs" in value["clinical"]["covariate_sets"]["mrs"]
    assert clinical.preflight_clinical(Path(__file__).resolve().parents[1], value)


def test_id_contract():
    assert clinical.normalize_ids(pd.Series(["0001","a002"]), {"kind":"string"}).tolist() == ["0001","a002"]
    assert clinical.normalize_ids(pd.Series([1.,20.]), {"kind":"numeric","width":4}).tolist() == ["0001","0020"]
    with pytest.raises(ValueError):
        clinical.normalize_ids(pd.Series([1.5]), {"kind":"numeric","width":4})
    with pytest.raises(ValueError):
        clinical.normalize_ids(pd.Series([1.]), {"kind":"numeric"})
    with pytest.raises(ValueError):
        clinical.normalize_ids(pd.Series([" A"]), {"kind":"string"})


def test_sas_special_missing(monkeypatch):
    fake = SimpleNamespace(read_sas7bdat=lambda *args, **kwargs: (pd.DataFrame({"score":[1,"A","_"]}), SimpleNamespace(missing_user_values={"score":["A","_"]}, column_names_to_labels={"score":"分数"}, file_label="test")))
    monkeypatch.setitem(sys.modules, "pyreadstat", fake)
    frame, metadata = clinical.read_sas(Path("unused.sas7bdat"))
    assert frame.score.isna().tolist() == [False,True,True]
    assert metadata["special_missing"]["score"] == ["A","_"]


def test_join_append_and_duplicates(tmp_path, config, monkeypatch):
    make_root(tmp_path)
    for name in ["data","batch","followup"]:
        (tmp_path / f"clinical/discovery/{name}.sas7bdat").write_bytes(b"fixture")
    records = {"data":pd.DataFrame({"id":["001"],"AGE":[50]}), "batch":pd.DataFrame({"id":["002"],"AGE":[60]}), "followup":pd.DataFrame({"id":["001","002","003"],"MRS":[1,2,3]})}
    monkeypatch.setattr(clinical, "read_sas", lambda path:(records[path.stem].copy(), {}))
    config["clinical"]["cohorts"]["discovery"]["id"] = {"field":"id","kind":"string"}
    frame, _ = clinical.load_sources(tmp_path, config, "discovery", tmp_path / "out")
    assert frame.subject_id.tolist() == ["001","002","003"]
    assert "AGE" in frame and "MRS" in frame
    records["batch"]["id"] = "001"
    with pytest.raises(ValueError, match="重复|冲突"):
        clinical.load_sources(tmp_path, config, "discovery", tmp_path / "bad")


def test_unmeasurable_and_complete_sum(tmp_path, config):
    frame = raw_frame(config).rename(columns=lambda name: "baseline::"+name)
    frame["subject_id"] = ["0001","a002","0003"]
    frame.loc[0,"baseline::A_NIH5A"] = 9
    mapped = clinical.map_variables(frame, config, "discovery", tmp_path)
    mapped["lesion_volume_mm3"] = [10,20,30]
    derived = clinical.derive_variables(mapped, config["clinical"]["derived"])
    assert np.isnan(derived.loc[0,"nihss_motor_left"])
    assert derived.loc[1,"nihss_motor_left"] == 0
    assert derived.loc[0,"log_lesion_volume"] == np.log1p(10)
    frame.loc[0,"baseline::A_NIH7"] = 9
    with pytest.raises(ValueError, match="编码或范围"):
        clinical.map_variables(frame, config, "discovery", tmp_path)


def test_ambiguous_field_is_not_first_column(tmp_path, config):
    frame = raw_frame(config).rename(columns=lambda name:"baseline::"+name)
    frame["subject_id"] = ["0001","a002","0003"]
    frame["other::AGE"] = 99
    with pytest.raises(ValueError, match="明确 source"):
        clinical.map_variables(frame, config, "discovery", tmp_path)
    config["clinical"]["cohorts"]["discovery"]["variable_overrides"]["age"] = {"source":"baseline"}
    assert clinical.map_variables(frame, config, "discovery", tmp_path).age.iloc[0] == 50


def test_death_conflict_and_structural_missing(tmp_path, config):
    config["statistics"]["families"]["mrs_single"]["enabled"] = True
    config["statistics"]["families"]["mrs_single"]["outcomes"] = ["mrs_3m","mrs_6m"]
    config["statistics"]["families"]["moca"]["enabled"] = True
    raw = pd.DataFrame({"subject_id":["001","002"]})
    mapped = pd.DataFrame({"subject_id":["001","002"],"mrs_3m":[6,1],"mrs_6m":[np.nan,1],"moca_3m":[np.nan,25],"moca_6m":[np.nan,25]})
    checked = clinical.check_death(raw, mapped, config, tmp_path)
    assert checked.moca_3m_structural_missing.tolist() == [True,False]
    assert np.isnan(checked.moca_3m.iloc[0])
    mapped.loc[0,"mrs_6m"] = 1
    with pytest.raises(ValueError, match="死亡状态"):
        clinical.check_death(raw, mapped, config, tmp_path)
    with pytest.raises(ValueError, match="重叠"):
        clinical.assert_disjoint({"001","002"}, {"002","003"})


def test_prepare_preserves_id_and_observed_outcomes(tmp_path, config, monkeypatch):
    make_root(tmp_path)
    original = raw_frame(config)
    original.loc[1,"A_NIH1A"] = np.nan
    monkeypatch.setattr(clinical, "read_sas", lambda path:(original.copy(), {}))
    features_dir = tmp_path / "features"
    features_dir.mkdir()
    frame = pd.DataFrame({"subject_id":["0001","a002","0003"], **{f"{channel}_{mode:03d}":[0.01,0.02,0.03] for channel in ["g","w"] for mode in range(1,130)}})
    frame.to_csv(features_dir / "mode_features.tsv.gz", sep="\t", index=False)
    pd.DataFrame({"subject_id":frame.subject_id,"lesion_volume_mm3":[10,20,30]}).to_csv(features_dir / "summary.tsv", sep="\t", index=False)
    manifest = clinical.prepare_clinical(tmp_path, config, features_dir, tmp_path / "out")
    data = pd.read_csv(tmp_path / "out/analysis.tsv.gz", sep="\t", dtype={"subject_id":str})
    assert data.subject_id.tolist() == ["0001","a002","0003"]
    assert data.age.isna().sum() == 1 and data.nihss_1a.isna().sum() == 1
    assert manifest["n"] == 3 and len(manifest["feature_columns"]) == 258
    assert data.log_lesion_volume.iloc[0] == pytest.approx(np.log1p(10))
    assert original.AGE.isna().sum() == 1
