#!/usr/bin/env python
"""总控的独立 Python 阶段工作进程；使临床导入/空间注释也能安全立即中断。"""
from pathlib import Path
import argparse
from functions.storage import read_json,save_json


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job')
    args=parser.parse_args()
    job=read_json(Path(args.job))
    root=Path(__file__).resolve().parent
    target=Path(job['output_dir'])
    config,upstream=job['config'],job['upstream']
    if job['stage']=='clinical':
        from functions.clinical import prepare_clinical,preflight_clinical,assert_disjoint
        from functions.cohort import configured_cohorts
        cohorts={}
        for cohort in configured_cohorts(config):
            issues=preflight_clinical(root,config,cohort)
            if issues:
                save_json(target/f'{cohort}_missing.json',{'issues':issues})
                if cohort=='discovery': raise ValueError('; '.join(issues))
                continue
            feature_dir = Path(upstream['lesion']) / cohort
            cohorts[cohort]=prepare_clinical(root,config,feature_dir,target/cohort,cohort)
        if 'validation' in cohorts and 'discovery' in cohorts:
            assert_disjoint(set(cohorts['discovery']['cohort_subject_ids']),set(cohorts['validation']['cohort_subject_ids']))
            # 记录两队列可共同验证的分析变量；R 阶段据此过滤验证结局。
            discovery_names = set(cohorts['discovery'].get('available_variables', []))
            validation_names = set(cohorts['validation'].get('available_variables', []))
            save_json(target/'common_variables.json',
                      {'variables': sorted(discovery_names & validation_names),
                       'discovery_only': sorted(discovery_names-validation_names),
                       'validation_only': sorted(validation_names-discovery_names)})
        save_json(target/'cohorts.json', cohorts)
    elif job['stage']=='annotation':
        from functions.spatial_export import prepare_combination_maps
        from functions.annotation import run_annotation
        maps=target/'maps'
        prepare_combination_maps(root,config,upstream,Path(job['clinical_dir']),maps)
        run_annotation(root,config,maps,target)
    else: raise ValueError('不支持的 Python 阶段')


if __name__=='__main__': main()
