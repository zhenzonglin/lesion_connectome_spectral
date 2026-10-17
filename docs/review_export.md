# 只导出汇总结果的审阅包

`export_review.py` 用于将当前运行的统计汇总和健康参考交给报告审阅者。它不运行分析、不修改配置、不写检查点、不访问原始 SAS、逐人病灶或患者分析宽表。脚本可放在轻量补丁目录中独立执行，依赖工作站现有环境中的 NumPy、PyYAML；无需安装补丁或停止总控。

## 工作站命令

在项目根目录运行，先检查预览，再使用同一个新文件名导出：

```bash
cd /data/usersdir/linzhenzong/lesion_connectome_spectral || exit 1
LCS_EXPORT=/data/usersdir/linzhenzong/lesion_connectome_spectral_hotfix_20261014/export_review.py
LCS_REVIEW_ZIP="$PWD/review_export_$(date +%Y%m%d_%H%M%S).zip"

./env/bin/python "$LCS_EXPORT" \
  --project "$PWD" --config config.workstation.yaml \
  --run main_20261014 --output "$LCS_REVIEW_ZIP" --dry-run

./env/bin/python "$LCS_EXPORT" \
  --project "$PWD" --config config.workstation.yaml \
  --run main_20261014 --output "$LCS_REVIEW_ZIP"
```

`--output` 必须指定，父目录必须已经存在。ZIP 独占新建，已有同名文件不会覆盖。导出过程中异常或源文件发生变化时，删除本次创建的未完成 ZIP，保留其他文件。`--dry-run` 验证同样的数据合同，仅打印文件数、字节数、列名、行数和缺失项目，不创建文件。实际导出成功时只打印ZIP绝对位置、字节数、文件数、reference_id、缺失项数量及枚举清理计数；完整清单保存在ZIP内。

## 来源和完整性

1. 只读取所指定 run 的 `stages.json` 当前指针。没有 archive、历史运行或测试结果回退。
2. 只导出状态为 `completed` 的 clinical、imputation、association、patterns、validation 阶段。失败／运行中／过期阶段记入 `missing_artifacts`，不将半成品视为完成结果。已完成阶段中的 `not_estimable` 模型行继续保留。
3. 阶段路径必须位于本 run 的 `stages/<stage>/<signature前16位>`；相同名称的其他目录、符号链接越界、完成签名不符均拒绝。
4. 每个选中文件必须在该阶段 `complete.json` 中，SHA256 必须一致。元数据从该阶段冻结的 `job.json` 提取；当前配置仅用于定位 results 根目录。
5. 导出前后再次检查所有实际读取的文件（包括状态指针）未变化。总控在导出时更新状态也可能触发安全中止；稍后重试即可，不需要重新计算。
6. 这是选定汇总文件的完整性检查，不会重新验证未导出的逐人文件，不会调用模型重拟合，也不能据此宣布所有模型或临床解释成立。

ZIP 的 `manifest.json` 记录阶段签名、相对阶段引用、导出文件 SHA256／行数／列名、缺失列表及未知枚举数量。run 名只保存摘要，不包含绝对目录。manifest 不自签名。

## 导出白名单

所有 CSV 均逐列重建。额外列、任意配置文本、自由文本原因和警告都不导出。当前版本仅支持项目已知 NIHSS 分项／合并项、七个单时点 mRS、纵向 mRS、三个 MoCA 结局和固定 G/W 特征。未知关键标识拒绝导出；未知状态枚举替换为 `unknown` 并计数，原行仍保留。数值必须可严格解析，概率在 `[0,1]`，模态在1–129；非法值不会被强制改成0。

| 数据 | 导出字段 |
|---|---|
| `association/associations.csv` | outcome、family、mode、analysis、status、n、observations、p_joint、mi_status、q、planned_family_tests；G/W各自estimate、standard_error、lower、upper、p |
| `association/jmax_sensitivity.csv` | 同上，但仅Jmax参数；未启用则记录缺失 |
| `association/threshold_sensitivity.csv` | outcome、family、mode、analysis、threshold、status、p_joint、term、estimate、standard_error、lower、upper、p |
| `association/time_interaction.csv` | outcome、family、mode、analysis、status、p_joint、q、planned_family_tests |
| `association/model_diagnostics_summary.csv` | 按outcome/family/mode/analysis汇总拟合数、PO状态计数、稀疏类别拟合数、最大条件数／梯度／天花板比例、最小诊断p；不输出每条警告文本 |
| `clinical/<cohort>/audit/sample_flow.csv` | outcome、family、image_matched、excluded_outcome、observed_outcome、missing_outcome、covariate_missing |
| `imputation/diagnostics.csv` | outcome、status、m、iterations、maximum_psrf、events |
| `imputation/imputation_summary.csv` | outcome、status、n |
| `patterns/patterns.csv` | target、status、lambda、lambda_boundary |
| `patterns/<target>/coefficients.csv` | feature、mean、scale、coefficient、raw_coefficient、selection_frequency、sign_agreement_with_ensemble、coefficient_sd |
| `patterns/<target>/cross_validation.csv` | lambda、mean_log_loss、standard_error、selected、converged_all_folds |
| `validation/validation.csv` | pattern、outcome、family、status、n、p_joint、estimate、lower、upper、training_diagnostics_passed、coefficient_direction_stable、q、planned_family_tests |
| `validation/pattern_eligibility.csv` | 上述字段加annotation_eligible与白名单枚举annotation_reason |
| `validation/planned_tests.csv` | pattern、outcome、family |
| `outcomes.json` | 按冻结阶段记录结局名、类型、家族、等级、协变量、纵向访视、规范方向与系数量纲；MI/FDR/模型/模式参数及死亡策略的有限数值和枚举；MI方法以及四个NIHSS合并项／体积log变换的本轮实际派生操作与规范输入。不包括SAS字段、自由文本标题、源路径、ID规则或配置全文 |

诊断表中的最小p只表示各拟合的警告范围，不是合并检验。validation现有CSV没有SE，本工具不从CI反推。阈值表可能有无穷自由度，导出合同不包含df，仅保留原有estimate、SE、CI和p。

## 健康参考

只读取 `reference/manifest.json` 指定并校验的 `modes.npz`、`sc.npz`、`nodes.tsv`，然后重建：

- `reference/healthy_modes.npz`：`mode_ids[129]`、`values[129]`、`vectors[135,129]`、`active[135]`、`strength[135]`、`weights[135,135]`。不复制未知NPZ成员，禁用pickle。SC权重和节点强度保留，后续可以解释W的健康边支撑；无需导出TRK/H5索引或患者断连数据。
- `reference/nodes.csv`：固定公开图谱的index、label、active、group、atlas_name、coverage_type、network_id及MNI质心坐标。严格验证135个已知标签顺序，防止患者ID混入允许的字符串列。
- `reference/spectrum.csv`：依据已校验健康数组重新汇总129个非零模态的特征值、有效节点数和皮层／皮层下／小脑／脑干平方载荷；不包含零模态。

核查135节点、130活动节点、129非零模态、1375边，以及SC对称性、非负权重、节点强度、谱顺序、正交性和特征方程。统计job的reference_id必须与当前健康参考一致。此处仅做既有健康参考的验证和汇总，不重新求特征向量。

## 明确排除

不打包：原始SAS、NIfTI、患者NPZ、`analysis.tsv.gz`、逐人特征／覆盖率／状态／排除表、`discovery_scores.csv`、`folds.csv`、`participants.csv`、任何RDS及MI填充值／链对象、含ID的审计、临床manifest、配置全文、job全文、日志、SQLite、HTML和既有患者图片。

尤其不能直接分享旧 `report.html`：其临床部分可能嵌入死亡冲突审计的患者ID。后续HTML应使用本工具生成的汇总数据重新制作。

## 解释限制

本包保留完整汇总行，支持核对方向、效应、SE、模型警告、样本量及全家族FDR分母；不会按显著性删行。G/W联合p没有单一方向，需看各自系数。原始G/W单位的效应未经标准化；连续MoCA系数与有序模型log odds不能直接混为同一量纲。`beta/SE`仍受样本量影响。

没有逐人数据，因此不能从本包执行患者bootstrap、患者亚型聚类或重新拟合模型；也没有跨结局／跨模态估计协方差。截图和本包中的候选组应明确标记为探索性，不能宣称已完成独立验证。需要稳定性分析时，应在工作站本地运行，另行导出聚合结果。

## 验证

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./env/bin/python -B -m pytest tests/test_export_review.py -q -p no:cacheprovider
```

测试只用临时目录中的合成汇总表与合成健康图，覆盖列注入、允许列中的秘密字符串、未知枚举保留行、非法数字、路径／符号链接越界、源哈希不符、源修改中止、禁止覆盖、未完成ZIP清理、dry-run只读、健康参考合同和多阶段聚合输出。合成测试不代表工作站真实数据已导出或审阅。
