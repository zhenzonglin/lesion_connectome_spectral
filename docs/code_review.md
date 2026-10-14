# 正式代码核查清单：按输入、处理、输出顺序阅读

这份清单审核代码和数据合同，不是患者病灶图片审核表。真实 SAS 尚未配置；临床模块通过合成测试只能证明程序合同，不能产生临床发现。运行状态和验收日志应与本清单同时阅读。

## 0. 先确定本次运行使用了什么

入口是 `run.py`，主配置是 `config.yaml`。`functions/configuration.py::load_config` 用 PyYAML 安全读取配置；`validate_contract` 直接拒绝不受支持的固定方法修改。代码不执行配置中的 Python/R 表达式。

可修改的内容包括数据路径、文件名/ID规则、来源字段和编码、启用的统计家族、协变量字典、插补数量与迭代、模型诊断阈值、交叉验证与惩罚强度、空间零模型数量及运行并行度。固定内容包括 LQT135 资源合同、G/W 定义、观察结局纳入规则、训练折隔离、D1/Rubin、BH 家族边界、原子提交以及不生成患者病灶质控图。固定内容仍写入配置，目的是可审核；不能把它们改为其他值而继续假称同一种方法。

每次正式运行目录是 `paths.results/<run>/`，默认 `results/<run>/`。阶段产物位于 `stages/<stage>/<指纹前16位>/`。指纹包含对应参数、输入内容、核心代码和上游完成记录。修改统计设置应使相关下游阶段失效，不让已完成病灶查询失效。历史与合成测试使用 `archive/...`。

| 先看哪里 | 核查内容 | 已固定 / 需根据数据修改 |
|---|---|---|
| `config.yaml::paths,subjects,imaging` | 输入目录、命名模式、实际MNI空间、1 mm、体积阈值 | 空间声明需用户真实确认；模板和135/130/129/1375校验值固定 |
| `config.yaml::clinical` | 真实 SAS、ID列、merge操作、变量、死亡和MoCA规则 | 候选字段必须人工核对；不能直接将 confirmed 改为 true 代替核查 |
| `config.yaml::statistics,patterns,validation,annotation` | 四个FDR家族、20目标、模型、MI、CV、验证与空间设置 | 科学方法合同固定；数值参数变更会进入新指纹 |
| `<run>/stages.json`、阶段 `parameters.json/complete.json` | 当前阶段状态、参数、实际产物SHA256 | completed表示阶段计算提交，不代表所有检验显著或可估计 |

## 1. 原始病灶与 ID

| 函数位置 | 输入 | 包与具体处理 | 输出位置 | 审核重点 |
|---|---|---|---|---|
| `functions/cohort.py::discover` | 病灶目录、`subjects.pattern` 或 participants TSV | pathlib/正则表达式；恰好一个 `*` 内全文作ID；查重复，计算SHA256 | lesion阶段 `manifest.json` 的被试清单 | 保留大小写和前导零；不去ses；一个ID只能一份病灶 |
| `functions/patient.py::process_subject` 空文件分支 | 文件字节数、只读原图 | 标准库先识别0字节文件；有效全零掩膜在空间读取后识别 | `subjects/<id>/record.json`，汇总 `excluded.tsv` | 状态必须是 excluded_empty_lesion；无result.npz、无伪造零损伤特征 |
| `functions/space.py::map_patient` | 病灶NIfTI、冻结参考template | NiBabel/NumPy；三维、二值、有限、1mm与仿射检查；必要时最近邻重采样 | 函数内存中的参考网格掩膜；record中记录前后体积与处理 | 数组下标不能代替物理坐标；未确认空间的文件不可凭1mm直接接受 |

病灶体积从原始有效掩膜的阳性体素数乘以仿射矩阵3×3部分行列式绝对值获得，单位 mm³。临床协变量使用这个原始体积；重采样后体积仅供空间一致性检查。体积相对变化大于 `imaging.maximum_volume_change` 阻止该例通过。代码不清理小灶或对侧体素，不保存默认患者质控PNG。

## 2. 冻结健康参考和脑区覆盖

| 函数位置 | 输入 | 包与处理 | 输出 | 固定 / 可改 |
|---|---|---|---|---|
| `functions/reference.py::verify_reference` | `reference/manifest.json` 及文件 | hashlib；逐文件核对SHA256 | 通过或明确异常 | 冻结参考不可用普通参数更换 |
| `load_reference`、`numerical_checks` | SC、模态、脑区表、模板、索引 | NumPy/SciPy/NiBabel；核对图对称、无自环、135全节点、130活动节点、129模态、1375边；谱残差与正交误差 | 内存reference、健康检查记录 | 零度排除清单保留；不得补造连接或零模态 |
| `functions/space.py::roi_damage` | 统一网格病灶、图谱、nodes表 | NumPy按标签统计覆盖：病灶∩脑区体素数/脑区总体素数 | 135维roi_coverage、命中体素数和解剖统计 | 覆盖率独立输出，不乘入SC；脑干分区为混合组织 |

健康 SC 来自同一 LQT/HCP842 流线底座，连续权重为流线计数除以两端脑区体积之和。健康基底使用对称归一化拉普拉斯。正式运行复用已经冻结并校验的结果，不逐例重构健康图，也不逐例重归一化。依据：[LQT作者代码](https://github.com/jcgriffis/lesion_quantification_toolkit)、[HCP842](https://doi.org/10.1016/j.neuroimage.2018.05.027)、[Preti 2019](https://doi.org/10.1038/s41467-019-12765-7)。

## 3. 病灶到流线中断

`functions/index.py::load_index` 读取压缩HDF5的体素—流线索引及固定端点分类。`query_ids` 用Numba查询所有病灶体素的流线编号并取并集。`interrupted_counts` 将有效不同脑区端点对汇总到1375条健康无向边。每条流线只计一次；同脑区、端点未匹配等类别保留统计，不冒充有效断连边。

运行时不需要原始TRK/FIB或DSI Studio。完整线段相交规则在一次性索引建立中已经编码；`append_segment/collect_memberships` 为构建和验收所用函数，正式患者查询不退回稀疏点采样。中断计数必须逐边满足0至健康计数上界。

结果保存在每例 `result.npz::lost_counts`，不默认保存中断TRK、1375×1375矩阵或重复完整SC。计数相减后除以冻结体积和可以按需恢复患者SC。

## 4. G、W、Jmax 与空间守恒

| 函数 | 输入 | 运算与含义 | 默认保存 |
|---|---|---|---|
| `functions/spectral.py::tissue_loss` | 135维覆盖率rho、135×129健康u | `rho[:,None]*u**2`后按节点求和；得到各模态平方载荷所在组织受损比例G | 每例129维g |
| `fixed_support_loss` | 1375边中断计数、体积、健康节点强度、u、lambda | 每条边的权重损失乘以健康归一化载荷差平方，边求和后除以对应lambda | 每例129维w，沿用原support_loss数值 |
| `process_subject` | g、w | 若 `save_jmax` 开启，逐模态取最大值；不是相加 | 每例129维jmax，可关闭保存 |
| `functions/cohort.py::export_tables` | 所有已提交非空患者 | pandas只汇总成功病例；错误/排除另列 | mode_features.tsv.gz、roi_coverage.tsv.gz、summary.tsv、status.tsv、excluded.tsv |

G/W均须在数值容差内处于[0,1]。按需请求贡献时，G节点贡献是覆盖率乘平方载荷；W把每条边的损失能量一半分给两个端点。节点贡献求和分别等于G/W。模态符号翻转不改变这些结果。未参加分解的五个节点在空间输出中为NaN，不解释为零损伤。G和W是不同损伤通道；不能把整个脑干都叫纯灰质，也不能将完整病灶中断的W称为仅白质掩膜分析。

G/W/Jmax属于本项目的损伤指标扩展，Preti文献提供健康谐波定义，不能引用该文献声称它已验证本项目的临床损伤指标。

## 5. SAS 到固定临床宽表

详细字段和代码逐步表见 `docs/clinical.md`。顺序为 `preflight_clinical → load_sources/normalize_ids/read_sas → map_variables → check_death → derive_variables → prepare_clinical`。主要包是pyreadstat、pandas、NumPy。

SAS导入保留标签和特殊缺失元信息；不根据字段顺序合并，不默认取重复第一条。每个来源使用明确base/append/join操作。发现与验证的原始队列成员做重叠检查，不能只检查影像成功者。单独的NIHSS `9`不可测编码按字段处理。死亡冲突报错，死后MoCA为结构性缺失，均不填补结局。

输出在临床阶段的 `discovery/analysis.tsv.gz` 和 `validation/analysis.tsv.gz`，列包括subject_id、258个g/w特征、原始体积和内部临床变量。`audit/`保存缺失、来源冲突、未匹配与样本流向。一个结局缺失不删除该人在其他结局的记录；筛选发生在下一步每个结局内。

## 6. 推断用多重插补

| R函数位置 | 输入 | 处理 | 产物及审核项 |
|---|---|---|---|
| `statistics/common.R::read_input,prepare_types` | analysis.tsv.gz、类型字典 | 保持ID文本；连续转数值；分类按明确levels编码 | 内存模型数据；不修改原始SAS |
| `statistics/imputation.R::endpoint_data` | 一项结局、协变量、完整G/W | 仅该结局已观察且影像完整者；必要协变量整列缺失/常量明确报错 | 固定患者名单，129模态复用 |
| `imputation_frame,mi_method` | 固定名单 | 推断MI含观察结局、G/W和其他协变量；纵向使用宽表观察状态辅助；只有协变量where允许填补 | methods.csv、predictor_matrix_requested.csv |
| `build_imputations` | 配置m、初始/递增/最大迭代、种子 | mice::mice/mice.mids；保存批次检查点；pmm/logreg/polyreg/polr | compact imputations.rds、participants.csv、actual predictor matrix、logged_events.csv |
| `rank_split_psrf,mi_diagnostics` | 各链均值/方差轨迹 | 秩归一化/分链PSRF；不完整诊断不算通过；超阈值递增迭代到上限 | chain_traces.rds、chain_diagnostics.csv、diagnostics.json |
| `complete_covariates` | 原始表、紧凑填充值、第j套MI | 只在记录的缺失协变量位置替换 | 单套内存数据，不保存30份完整G/W表 |

默认m=30、初始20、递增20、最多100、PSRF警报1.1。达到上限仍有问题保留警报，不退回单次插补。结局、影像和体积从不填补。PSRF是诊断工具，不能单独证明插补模型有效。纵向宽表FCS与有序随机效应模型的相容性是近似的，报告必须保留此限制。依据：[MICE官方文档](https://amices.org/mice/reference/mice.html)。

## 7. 逐模态模型与四个FDR家族校正

`statistics/inference.R::fit_model` 对每份MI分别拟合：NIHSS/mRS用ordinal::clm；MoCA用lm后sandwich::vcovHC(HC3)；纵向mRS用ordinal::clmm，时间分类、患者随机截距。每个模态联合放入G和W，协变量按对应结局集加入。prestroke_mrs与教育用分类编码。

代码检查有限估计、信息矩阵、梯度/优化收敛和稀疏等级。有序横断面使用nominal_test诊断比例优势；每套MI诊断分别保留，不平均这些p。比例优势警报触发预设阈值敏感性；纵向不宣称存在原生nominal_test，其敏感性用glmer保留患者随机效应。MoCA保存残差分布、异方差和天花板比例。

`pool_parameters/pool_fits` 合并G/W两参数的D1联合检验，单个系数用Rubin规则与自由度计算区间；禁止平均p、SE或似然比。`analyze_terms` 汇总假设和MI警报；不可估计项目记录原因。

`association_stage` 以结局/模态写 `mode_*.rds` 检查点和汇总 `associations.csv`。个体系数、协方差与敏感性保存在检查点中。`statistics/common.R::family_adjust` 按四个FDR家族BH，计划检验数包含失败项，失败p/q保留缺失而不是当作不显著。

| 家族 | 预设检验范围 | 默认模型 |
|---|---|---|
| nihss | 19结局×129模态 | 有序logistic |
| mrs_single | 7月份×129模态 | 有序logistic |
| mrs_longitudinal | 整体纵向×129模态，默认关闭 | 有序随机截距模型 |
| moca | 3/6/12月共用129模态家族 | 线性+HC3 |

依据：[D1](https://amices.org/mice/reference/D1.html)、[mitml多参数约束](https://search.r-project.org/CRAN/refmans/mitml/html/testConstraints.html)、[ordinal](https://search.r-project.org/CRAN/refmans/ordinal/html/ordinal-package.html)。

## 8. 20个组合学习与防泄漏

`statistics/patterns.R::fit_pattern` 分别学习mRS3月和NIHSS19结局。全部258影像特征进入候选，不按单模态显著性筛选。`stratified_folds` 在患者层面分折；`image_scaler` 只由训练行估计均值、标准差、常量特征。`pattern_matrix` 将标准化影像与不惩罚协变量连接。

训练用MI调用 `build_imputations(inference=FALSE, ignored=验证行)`，不把目标结局放入插补器，验证行不参与插补参数估计。`fit_net` 调用ordinalNet累计logit弹性网，影像penaltyFactors=1、协变量=0；alpha、lambda网格和收敛参数均记录。

每个lambda平均30套MI的验证概率，再计算每折log loss；标准误来自患者折，不把MI数当独立样本。一标准误规则选择符合阈值的最大lambda；网格边界被标记，不自动试更多网格后挑结果。全发现集重拟合后在同一尺度平均影像系数，得到固定影像得分；协变量只用于学习调整，不写入影像得分。惩罚后系数不套普通Rubin显著性检验。

输出为每个目标的 `folds.csv`、`cross_validation.csv`、`coefficients.csv`、`imputation_coefficients.rds`、`model.json`、`discovery_scores.csv` 和恢复检查点；阶段汇总为 `patterns.csv`。model冻结特征顺序、标准化、系数、参考版本、目标方向与训练签名。依据：[ordinalNet接口](https://search.r-project.org/CRAN/refmans/ordinalNet/html/ordinalNet.html)。

## 9. 外部验证

`frozen_score` 读取模型与验证G/W，仅执行冻结的影像线性组合，不重新估计权重、尺度或惩罚参数。`validation_stage` 先写 `planned_tests.csv`，对各模式指定结局做调整关联；主要mRS组合额外关联七时点mRS、19项NIHSS、三时点MoCA。验证集协变量独立MI，可用已观察结局辅助推断，因此这些结果是关联验证，不是未知结局的部署预测性能。

`validation.csv` 按外部验证预设家族校正；`pattern_eligibility.csv` 标记是否具备后续机制注释资格。未复现、方向不一致或严重假设警报不能自动升级。发现集和验证集分别校正，不能合并后挑显著模式。

## 10. 三类空间图与注释

| 函数位置 | 输入、包、处理 | 输出 | 解释限制 |
|---|---|---|---|
| `functions/annotation.py::export_healthy_modes` | 冻结u；NumPy、NiBabel、Matplotlib | 129模态各signed/squared目录，regions.tsv与map.png | 健康背景，不是患者功能激活 |
| `loss_contributions` | 每例计数/覆盖；调用固定G/W公式 | 按需135×129内存贡献 | 不默认保存大矩阵 |
| `combination_maps` | 单例贡献、已冻结258系数/均值/尺度 | 分开的G/W节点图、score、centering_constant、守恒误差 | 原始贡献总和减中心化常数=标准化得分；绝不以sum(a*u)替代 |
| `mean_combination_maps` | 指定队列紧凑NPZ | 在线平均两通道贡献，不缓存患者×节点×模态 | 平均图受指定队列病灶分布影响，应标明队列 |
| `load_annotations` | 资源manifest与图谱/节点SHA | 校验Schaefer对应和资源版本 | 不能只凭100行直接拼表 |
| `moran_surrogates,correlation_test` | 有效皮层图、半球/坐标、固定术语/受体矩阵 | 空间零模型与Pearson相关检验 | 当前分半球Moran谱随机化；并非随机打乱脑区 |
| `run_annotation` | 预设maps_manifest、验证状态 | associations.tsv/png、null_diagnostics.json、summary.json | 两资源分别校正全部预设模式×G/W×术语；未验证对象不形成机制结论 |

正式推断首版限定正确匹配、完整覆盖的100个皮层节点。皮层下、小脑、脑干只报告覆盖，不能套用皮层空间零模型。Neurosynth/Hansen相关只用于解释已确定的模式，不反向训练临床权重。资源依据：[Hansen作者数据](https://github.com/netneurolab/hansen_receptors)、[空间零模型说明](https://github.com/netneurolab/neuromaps/blob/main/docs/user_guide/nulls.rst)。

## 11. 提交、停止、续跑与离线包

| 函数或入口 | 行为 | 要审核的证据 |
|---|---|---|
| `functions/runner.py::coordinate` | SQLite记录病例；ProcessPoolExecutor患者并行；完整文件检查后提交 | 一个成功患者仅紧凑NPZ与record；状态与SHA一致 |
| `functions/pipeline.py::run_pipeline` | 按阶段依赖和指纹执行；下游缺依赖明确阻止；完成清单核对产物 | stages.json、complete.json、上下游签名；续跑不会混入变更输入 |
| `run_child`、R `check_stop` | 只管理本次创建的子进程；正常停止在病例/统计检查点处理 | 停止后状态stopped；未提交部分可重算，已提交结果不覆盖 |
| `monitor.py` | 默认每30秒持续读取状态和经身份校验的进程，`--once`单次查看 | 各阶段、双队列病例计数和模型当前任务；识别前次运行记录，失效PID不能当运行中 |
| `stop.py`、`functions/control.py` | 核查PID、创建时间、启动token、路径；立即中断只限本任务 | 不应使用pkill或终止WSL等广泛动作 |
| `package.py`、`setup.sh` | 代码、环境、冻结参考和注释资源离线移植 | 文件清单和SHA；默认排除SAS、病灶、结果、archive和缓存 |

原始数据、已验证历史结果的SHA必须保持不变。代码的通过状态、不可估计模型数量和真实数据是否就位分别报告，不把运行完成写成临床分析成功。

## 12. 人工代码审核记录

可在本表填写代码审核意见，不要求逐例病灶图审核。

| 核查项 | 状态（待核查/通过/需修改） | 审核人 | 备注 |
|---|---|---|---|
| 配置固定合同与实际可调参数一致 | 待核查 | | |
| ID、空病灶、体积和空间输入检查 | 待核查 | | |
| 健康参考、索引与G/W数值定义 | 待核查 | | |
| SAS合并、字段、死亡/MoCA规则 | 待核查 | | |
| 观察结局名单、MI和诊断 | 待核查 | | |
| D1/Rubin、PO、HC3、纵向模型 | 待核查 | | |
| 四个FDR家族预设检验清单和失败项 | 待核查 | | |
| 折内预处理、固定组合与独立验证 | 待核查 | | |
| 空间贡献守恒与注释零模型 | 待核查 | | |
| 阶段恢复、文件保护、离线包隐私排除 | 待核查 | | |
