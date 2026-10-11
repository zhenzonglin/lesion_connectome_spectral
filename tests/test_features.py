"""新增 G/J、输入命名与阶段独立性的工程验证。"""
from pathlib import Path
from copy import deepcopy
import numpy as np
import pytest
import nibabel as nib
from functions.reference import load_reference
from functions.spectral import tissue_loss,fixed_support_loss
from functions.cohort import discover
from functions.space import map_patient
from functions.pipeline import stage_parameters


def test_tissue_conservation_sign_and_exclusion():
    root = Path(__file__).resolve().parents[1]
    reference = load_reference(root,with_index=False)
    rho = np.linspace(0,1,135)
    local,spatial = tissue_loss(rho,reference,True)
    assert np.allclose(local,np.nansum(spatial,axis=0))
    changed = dict(reference,vectors=-reference['vectors'])
    assert np.array_equal(local,tissue_loss(rho,changed)[0])
    assert np.isnan(spatial[~reference['active']]).all()
    assert np.allclose(tissue_loss(np.ones(135),reference)[0],1)
    assert not tissue_loss(np.zeros(135),reference)[0].any()
    first,second = reference['edge_i'],reference['edge_j']
    white = fixed_support_loss(reference['counts'][first,second]*.2,reference)[0]
    joint = np.maximum(local,white)
    assert np.all(joint>=local) and np.all(joint>=white) and np.all(joint<=1)


def test_naming_preserves_exact_id(tmp_path):
    (tmp_path/'lesion').mkdir()
    for name in ['sub-0001_lesion.nii.gz','sub-Abses01_lesion.nii.gz']:
        (tmp_path/'lesion'/name).touch()
    config = {'paths':{'lesion':'lesion'},'subjects':{'pattern':'sub-*_lesion.nii.gz'},
              'space':'MNI152NLin6Asym','imaging':{'space_confirmed':True}}
    items = discover(tmp_path,'all',config)
    assert [item['subject_id'] for item in items] == ['0001','Abses01']
    with pytest.raises(ValueError,match='重复'): discover(tmp_path,'0001,0001',config)
    config['subjects']['pattern']='sub-**_lesion.nii.gz'
    with pytest.raises(ValueError,match='恰好'): discover(tmp_path,'all',config)


def test_resolution_rejected(tmp_path):
    image = nib.Nifti1Image(np.ones((3,3,3)),np.diag([2,2,2,1]))
    path=tmp_path/'wrong.nii.gz'
    nib.save(image,path)
    with pytest.raises(ValueError,match='1 mm'): map_patient(path,image)


def test_downstream_config_does_not_invalidate_imaging():
    config={'imaging':{'space':'MNI152NLin6Asym'},'subjects':{'pattern':'sub-*_lesion.nii.gz'},
            'clinical':{'variables':{}},'statistics':{'mi':{'m':30}}}
    changed=deepcopy(config)
    changed['statistics']['mi']['m']=40
    assert stage_parameters('lesion',config)==stage_parameters('lesion',changed)
    assert stage_parameters('imputation',config)!=stage_parameters('imputation',changed)
