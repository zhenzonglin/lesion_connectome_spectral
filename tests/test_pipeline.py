"""阶段缓存、原子提交与中断的独立工程测试。"""
import os
import time
import threading
import sys
from pathlib import Path
import pytest
from functions.control import run_lock,matching_process
from functions.pipeline import run_child,clinical_outcomes
from functions.storage import save_json,verify_files,digest
from functions.configuration import load_config


def test_repeated_lock_and_stale_identity(tmp_path):
    with run_lock(tmp_path):
        with pytest.raises(RuntimeError,match='重复'): 
            with run_lock(tmp_path): pass
    assert matching_process(tmp_path,{'pid':os.getpid(),'created':0,'root':str(tmp_path),'boot_id':'invalid'}) is None


def test_partial_output_is_not_completed(tmp_path):
    path=tmp_path/'result.txt'
    path.write_text('complete')
    expected={'result.txt':digest(path)}
    verify_files(tmp_path,expected)
    path.write_text('partial')
    with pytest.raises(ValueError): verify_files(tmp_path,expected)


@pytest.mark.parametrize('immediate',[False,True])
def test_stops_only_owned_child(tmp_path,immediate):
    # 子进程在批次边界看 stop 文件；立即停止由父进程终止本次进程组。
    script=tmp_path/'batch.py'
    script.write_text('import time\nfrom pathlib import Path\n'
                      'for i in range(100):\n time.sleep(.05)\n if Path(__file__).with_name("stop").exists(): break\n')
    state={'graceful':False,'immediate':False}
    request=tmp_path/'request.json'
    def stop():
        time.sleep(.15)
        save_json(request,{'token':'test','immediate':immediate})
    thread=threading.Thread(target=stop)
    thread.start()
    code=run_child([sys.executable,str(script)],tmp_path,state,request,'test',lambda:None)
    thread.join()
    assert state['graceful'] and code == (130 if immediate else 0)


def test_four_families_and_twenty_targets():
    root=Path(__file__).resolve().parents[1]
    config=load_config(root/'config.yaml')
    items=clinical_outcomes(config)
    assert len(config['statistics']['families'])==4
    assert sum(item['family']=='nihss' for item in items)==19
    assert sum(item['family']=='mrs_single' for item in items)==7
    assert len(config['patterns']['targets'])==20
    assert all('acute_nihss' not in item['covariates'] for item in items if item['family']=='nihss')
