"""将冻结组合映射到验证队列的平均节点贡献；不把健康载荷相加当成损伤图。"""
from pathlib import Path
import numpy as np
import pandas as pd
from functions.storage import read_json,save_json
from functions.reference import load_reference
from functions.annotation import combination_maps
from functions.plotting import mode_png


def prepare_combination_maps(root: Path, config: dict, upstream: dict, clinical_dir: Path,
                             output_dir: Path) -> dict:
    """先平均中断计数与覆盖率，再作线性贡献映射；等于逐例贡献的平均。"""
    output_dir.mkdir(parents=True,exist_ok=True)
    eligibility = pd.read_csv(Path(upstream['validation'])/'pattern_eligibility.csv')
    cohort = pd.read_csv(clinical_dir/'validation/analysis.tsv.gz',sep='\t',dtype={'subject_id':str})
    reference = load_reference(root,with_index=False)
    expected = [f'{channel}_{mode:03d}' for channel in ['g','w'] for mode in range(1,130)]
    maps,checks = [],[]
    for pattern in config['patterns']['targets']:
        own = eligibility.loc[eligibility['pattern']==pattern]
        eligible = len(own)==1 and str(own.iloc[0]['annotation_eligible']).lower()=='true'
        for channel in ['g','w']:
            maps.append({'pattern_id':pattern,'channel':channel,'map_file':f'{pattern}_{channel}.tsv','validated':eligible})
        if not eligible: continue
        model = read_json(Path(upstream['patterns'])/pattern/'model.json')
        if model['feature_columns'] != expected: raise ValueError('组合特征顺序与冻结健康模态不一致')
        if model['reference_id'] != read_json(root/'reference/manifest.json')['reference_id']:
            raise ValueError('组合与本项目健康参考版本不同')
        selected = cohort.loc[cohort[pattern].notna(),'subject_id']
        if not len(selected): raise ValueError('可注释模式没有对应结局的验证患者')
        lost_sum = np.zeros(len(reference['edge_i']),dtype=float)
        rho_sum = np.zeros(135,dtype=float)
        for subject in selected:
            with np.load(Path(upstream['lesion'])/'subjects'/subject/'result.npz') as data:
                lost_sum += data['lost_counts']
                rho_sum += data['roi_coverage']
        result = combination_maps(reference,lost_sum/len(selected),rho_sum/len(selected),
                                  model['coefficients'],model['feature_mean'],model['feature_scale'])
        checks.append({'pattern':pattern,'n':len(selected),'score_mean':result['score'],
                       'centering_constant':result['centering_constant'],
                       'conservation_error':result['conservation_error']})
        for channel in ['g','w']:
            pd.DataFrame({'index':np.arange(1,136),'value':result[channel]}).to_csv(
                output_dir/f'{pattern}_{channel}.tsv',sep='\t',index=False)
            mode_png(reference['template'],reference['atlas'],result[channel],
                     f'{pattern} | {channel.upper()} mean loss contribution',output_dir/f'{pattern}_{channel}.png',
                     note=f'external cohort observed target n={len(selected)}; signed raw-scale contribution',
                     dpi=config['annotation']['figure']['dpi'],
                     cmap_signed=config['annotation']['figure']['colormap_signed'],
                     cmap_loss=config['annotation']['figure']['colormap_loss'])
    manifest={'maps':maps,'object':'mean patient loss contribution in the observed-target external cohort',
              'checks':checks,'interpretation':'两通道节点和=平均影像得分+冻结中心化常数；不是空间因果效应'}
    save_json(output_dir/'maps_manifest.json',manifest)
    return manifest
