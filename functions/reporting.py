"""生成可离线打开的轻量报告；图片用相对路径，不嵌入重复大文件。"""

from pathlib import Path
import html
import pandas as pd
from functions.storage import read_json


def write_report(root: Path, directory: Path) -> None:
    progress = read_json(directory / "progress.json")
    summary = (
        pd.read_csv(directory / "summary.tsv", sep="\t")
        if (directory / "summary.tsv").stat().st_size > 1
        else pd.DataFrame()
    )
    from functions.plotting import cohort_spectrum

    spectrum = directory / "mode_features.tsv.gz"
    if len(summary):
        cohort_spectrum(spectrum, directory / "spectrum.png")
    document = """<!doctype html><html lang="zh"><meta charset="utf-8"><title>病灶结构谐波分析</title>
<style>body{font:16px system-ui;max-width:1100px;margin:32px auto;padding:0 20px;color:#183345;background:#f7fafb}h1,h2{color:#174b62}table{border-collapse:collapse;font-size:13px}td,th{padding:8px;border:1px solid #cbd8df}img{width:100%}section{background:white;padding:16px;margin:20px 0}code{background:#e6edf2;padding:3px}</style>"""
    document += f'<h1>LQT135 · 病灶结构谐波分析</h1><p>计算状态：<b>{html.escape(progress["status"])}</b>。成功 {progress.get("done",0)} 例，空病灶排除 {progress.get("excluded_empty_lesion",0)} 例，失败 {progress.get("failed",0)} 例。</p>'
    document += "<p>按用户要求，不进行患者病灶质控，不生成病灶叠加 PNG 或患者人工审核表。程序保留计算所需的文件格式、坐标及重采样一致性检查。空病灶记录在 excluded.tsv，不生成零值特征。</p>"
    document += "<p>135 个解剖分区，130 个健康有连接节点，129 个非零健康模态。灰质覆盖率独立输出；脑干为混合组织。未接入的 5 个节点不作为零损伤解释。</p>"
    document += "<p>网络指标表示冻结健康参考中的结构支撑损失，是本项目方法扩展；不直接等于神经活动、功能丧失或重组。患者形状分析保留在历史探索归档。</p>"
    document += "<p>HCP-842 为群体平均流线图谱；输入假定已完成 MNI 配准，本阶段没有进行患者解剖质量评定。未匹配端点和图谱外组织的覆盖不足会影响可估计范围。</p>"
    document += (
        summary.head(100).to_html(index=False, escape=True)
        if len(summary)
        else "<p>尚无通过计算的病例。</p>"
    )
    if len(summary):
        document += '<section><h2>健康模态支撑损失</h2><img src="spectrum.png" alt="mode support loss"></section>'
    document += "<p>汇总表显示前 100 例；完整清单见 summary.tsv、status.tsv。空病灶见 excluded.tsv；它们不进入模态或脑区覆盖率特征表。</p>"
    document += '<p>方法：<a href="https://doi.org/10.1016/j.nicl.2021.102639">LQT</a>；<a href="https://doi.org/10.1038/s41467-019-12765-7">Preti 2019</a>。临床关联和外部图谱注释尚未执行。</p></html>'
    (directory / "report.html").write_text(document, encoding="utf-8")
