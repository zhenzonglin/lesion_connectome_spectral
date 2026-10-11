"""汇集当前指纹的真实结果：样本流向、缺失、统计诊断、组合和注释；不替换缺失阶段。"""
from pathlib import Path
import html
import os
import json
import numpy as np
import pandas as pd


def table_section(path: Path, title: str, directory: Path, limit: int = 30) -> str:
    if not path.is_file(): return ''
    try:
        frame=pd.read_csv(path,sep='\t' if path.suffix=='.tsv' else ',')
    except pd.errors.EmptyDataError:
        return f'<section><h2>{html.escape(title)}</h2><p>没有可展示的记录；不生成零值替代。</p></section>'
    url=html.escape(os.path.relpath(path,directory),quote=True)
    return f'<section><h2>{html.escape(title)}</h2><p>共{len(frame)}行，展示前{limit}行。<a href="{url}">完整表</a></p>'+frame.head(limit).to_html(index=False,escape=True)+'</section>'


def report_sections(directory: Path, state: dict) -> str:
    """只读取已完成阶段的当前输出，不从archive借用测试或历史结果。"""
    parts=[]
    for name,record in state['stages'].items():
        if record.get('status')!='completed' or not record.get('path'): continue
        target=Path(record['path'])
        if name=='clinical':
            for cohort in ['discovery','validation']:
                for file,title in [('sample_flow','样本流向'),('missingness','缺失分布'),('death_conflicts','死亡冲突')]:
                    parts.append(table_section(target/cohort/'audit'/f'{file}.tsv',f'{cohort} · {title}',directory))
        elif name=='imputation':
            records=[]
            for path in sorted(target.glob('*/diagnostics.json')):
                value=json.loads(path.read_text())
                records.append({'outcome':path.parent.name, **{key:value.get(key) for key in
                    ['status','m','iterations','maximum_psrf']}})
            if records:
                parts.append('<section><h2>协变量插补诊断</h2>'+pd.DataFrame(records).to_html(index=False,escape=True)+'</section>')
        elif name=='association':
            path=target/'associations.csv'
            if path.exists():
                data=pd.read_csv(path)
                summary=data.groupby(['family','status'],dropna=False).size().reset_index(name='tests')
                parts.append('<section><h2>四家族及模型状态</h2>'+summary.to_html(index=False,escape=True)+'</section>')
                parts.append(table_section(path,'逐模态联合检验与效应',directory))
                diagnostic_path=target/'model_diagnostics.csv'
                if diagnostic_path.is_file():
                    diagnostics=pd.read_csv(diagnostic_path)
                    if {'outcome','po_status'}<=set(diagnostics):
                        diagnostic_summary=diagnostics.groupby(['outcome','po_status'],dropna=False).size().reset_index(name='fits')
                        parts.append('<section><h2>有序模型比例优势诊断</h2>'+diagnostic_summary.to_html(index=False,escape=True)+'</section>')
                    parts.append(table_section(diagnostic_path,'完整模型诊断（含 HC3 模型诊断）',directory,10))
                for filename,title in [('threshold_sensitivity.csv','阈值敏感性'),('jmax_sensitivity.csv','Jmax 次级分析'),('time_interaction.csv','纵向时间交互')]:
                    parts.append(table_section(target/filename,title,directory,10))
                if len(data) and {'outcome','mode','q'}<=set(data):
                    import matplotlib
                    matplotlib.use('Agg')
                    import matplotlib.pyplot as plt
                    matrix=data.pivot(index='outcome',columns='mode',values='q')
                    figure,axis=plt.subplots(figsize=(13,max(3,.3*len(matrix))),dpi=150)
                    values=-np.log10(matrix.to_numpy().clip(1e-12,1))
                    artist=axis.imshow(np.ma.masked_invalid(values),aspect='auto',cmap='viridis')
                    axis.set_yticks(range(len(matrix)),matrix.index,fontsize=8)
                    axis.set(xlabel='Healthy mode (ascending)',ylabel='Outcome',title='Joint G/W association | -log10 family q')
                    figure.colorbar(artist,ax=axis,label='-log10(q)')
                    figure.tight_layout()
                    figure.savefig(directory/'association.png')
                    plt.close(figure)
                    parts.append('<img src="association.png" alt="joint association">')
        elif name=='patterns':
            parts.append(table_section(target/'patterns.csv','组合学习状态',directory))
            for path in sorted(target.glob('*/coefficients.csv')):
                parts.append(table_section(path,f'{path.parent.name} · 冻结权重',directory,10))
        elif name=='validation':
            parts.append(table_section(target/'validation.csv','外部验证',directory))
            parts.append(table_section(target/'pattern_eligibility.csv','机制注释资格',directory))
        elif name=='annotation':
            parts.append(table_section(target/'associations.tsv','空间零模型注释',directory))
            image=target/'associations.png'
            if image.exists(): parts.append(f'<img src="{os.path.relpath(image,directory)}" alt="spatial annotation">')
    return ''.join(parts)
