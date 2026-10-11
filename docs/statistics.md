# 正式临床统计后端说明

本说明对应 `statistics/` 内的 R 程序。真实 SAS 尚未配置时，正式临床阶段应停在输入检查；`archive/statistics_*` 中的数据与结果均为工程测试，不构成临床发现。

## 1. 入口、输入与输出

启动方式为 `./env/bin/Rscript statistics/run.R job.json`。通常由 Python 总控生成 JSON，不需要手工维护两套分析配置。

JSON 必需字段：

| 字段 | 含义 |
|---|---|
| `stage` | `imputation`、`association`、`patterns` 或 `validation` |
| `input_csv` | Python 整合后的患者宽表，支持 CSV、TSV 及 `.gz` |
| `output_dir` | 当前统计阶段结果目录 |
| `variables` | 字段类型与分类变量的明确等级 |
| `outcomes` | 结局名称、类型、家族、等级、协变量；纵向结局另给访视列和月份 |
| `feature_columns` | 固定顺序的 258 个 G/W 列名 |
| `mode_pairs` | 129 个 `{mode,g,w}` 对应关系 |
| `config` | 主配置中的统计与组合学习设置 |

可选字段包括 `stop_file`、`upstream_imputation_dir`、`frozen_models_dir`、`reference_id` 和不可插补列 `immutable_columns`。JSON 不包含可执行 R 表达式。变量名只能是安全的字母、数字和下划线。

`subject_id` 始终为字符串，保留前导零、大小写和字面 ID `NA`。结果保存到 Python 指定的 `results/<run>/...` 阶段目录；测试指定 `archive`。入口核对输入校验值、配置和 R 源码，变化后拒绝向原目录混写。

## 2. 按代码顺序核查

| 文件与函数 | 输入 | 调用包与处理 | 输出与位置 |
|---|---|---|---|
| `common.R::read_input` | 标准临床宽表 | 基础 R 读入，检查 ID 与列名；ID 单独按字符串读取 | 内存宽表；不修改 SAS |
| `imputation.R::endpoint_data` | 宽表、一个结局的定义 | 固定该结局已观察且影像完整的名单；按配置转类型；检查必要协变量 | 端点固定患者表；缺失影像/体积不插补 |
| `imputation.R::build_imputations` | 固定名单、协变量与完整 G/W | `mice` 的 PMM/logreg/polyreg/polr，只对协变量开放 `where`；推断插补加入已观察结局 | `imputations.rds` 只含缺失位置的填充值；方法、预测矩阵、事件、链诊断与成员 CSV |
| `imputation.R::complete_covariates` | 原表、紧凑填充值、插补编号 | 在内存中还原一套协变量；已观察值和 G/W 保持不变 | 内存完整协变量表，不保存 30 份影像数据 |
| `inference.R::fit_model` | 一套协变量、G/W 或固定组合得分 | `ordinal::clm`、`ordinal::clmm` 或 `lm`+`sandwich::vcovHC(HC3)`；检查数值与假设 | 小型参数、协方差和诊断对象 |
| `inference.R::pool_parameters` | 每套参数估计及协方差 | `mitml::testConstraints(D1)` 联合检验；系数按 Rubin 合并；单参数用等价 Rubin t 检验 | 联合 p、系数、置信区间；不平均 p 或标准误 |
| `inference.R::association_stage` | 固定端点插补与 129 个模态 | 每模态拟合 G/W 联合模型；固定家族大小 BH 校正 | `associations.csv`、`model_diagnostics.csv`、每模态 RDS 检查点 |
| `patterns.R::fit_pattern` | 258 影像特征、一个有序目标及协变量 | `ordinalNet`，训练折内标准化和插补；协变量不惩罚；1SE 选择 | 每目标 `model.json`、全部权重、MI 系数、折结果及发现集得分 |
| `patterns.R::validation_stage` | 外部表、冻结模型 | 检查成员重叠；只应用原尺度权重；对固定得分进行协变量调整后的关联推断 | `planned_tests.csv`、`validation.csv`、`pattern_eligibility.csv` |

## 3. 缺失处理与插补诊断

主分析逐结局纳入已观察结局者。连续结局、mRS、NIHSS 原始分项或合并项本身均不插补。病灶体积及其对数变换来自影像，不能作为缺失临床协变量填补。必要协变量整列缺失或没有变异时，记录不可估计原因。

默认 30 套插补，初始 20 次迭代，按 20 次增加至最多 100 次。使用秩归一化、分链与折叠 PSRF，分别检查链均值和链方差，阈值 1.1。该诊断在基础 R 中实现，不需要额外安装 Stan。轨迹、实际预测矩阵、共线性事件均保留。PSRF 未定义或达到上限仍异常时，保留警报，不能当作通过。

插补链进行中的 `mi_state.rds` 是恢复检查点，含原始宽表一份；最终填充值原子提交后删除该临时链对象，只保留缺失协变量的填充值和小型链轨迹。

组合学习的插补与推断不同。每折首先**只用训练行**完成多重插补，不加入目标结局。然后对每套已完成的训练表，将其作为固定拟合行和供者，条件填补验证折协变量。训练侧填充值不会随验证侧标签、数值或缺失情况变化。验证侧每套训练插补配一个条件抽样，合起来仍有 30 套验证协变量，不能把单个条件抽样误解为整个分析只做一次插补。

## 4. 模型与假设

单时点 mRS、NIHSS 采用累积 Logistic 有序模型。程序将连续协变量做等价的数值重参数化以改善条件数；G/W 效应仍按原始损伤单位报告。分类协变量不当作等距连续数值。检查协方差正定性、条件数、优化器状态和信息矩阵尺度的梯度。

比例优势诊断使用 `ordinal::nominal_test`。调用前固定已展开公式和实际建模数据，防止 R 动态公式在更新模型时失效。每套插补的诊断 p 单独保存，不把诊断 p 平均成合并 p。发现假设警报后保留预设主模型，并输出所有阈值的二分类敏感性结果，不选择最有利阈值。二等级结局没有额外比例优势约束。

MoCA 使用线性模型与 HC3 协方差，再进行 MI 合并。输出残差分位数、异方差与天花板比例诊断，不将 MoCA 的缺失或死亡后未测评当作零分。

纵向 mRS 默认关闭。启用后采用分类时间与患者随机截距；允许部分访视缺失。协变量在患者宽表中插补，已观察等级和缺失状态作为完整辅助信息，再映射到真实访视。该 FCS 是近似方案，不宣称与混合模型完全相容。`clmm` 没有原生 `nominal_test`；程序输出保留患者随机效应的阈值敏感性模型。可选时间交互单独写入 `time_interaction.csv`，作为次级检验集合。

四个主家族为 NIHSS、七个单时点 mRS、纵向 mRS 和 MoCA（3、6、12 月共用一个 FDR 家族）。失败模型仍计入计划检验总数，其 p、q 为缺失而非 1 或 0。若启用 Jmax，程序另外拟合 `Jmax_k + 协变量`，保存 `jmax_sensitivity.csv`，按同样四个家族在独立的次级检验集合中校正；从不把 G/W/Jmax 同时纳入一个模型。

## 5. 组合、独立验证与注释资格

20 个组合均使用全部 258 个 G/W，不按逐模态显著性预筛选。5 折交叉验证；每折内独立估计影像均值、尺度与常量列。临床协变量的惩罚因子为 0，影像为 1。每个惩罚强度先平均 MI 验证概率，再计算患者折的对数损失；标准误来自折，而不是把插补次数当独立样本。

1SE 规则在最优损失的一标准误范围内选择最大的惩罚强度。网格边界和拟合失败必须报告。全发现集重拟合后，30 套影像系数取平均形成固定组合，属于 MI 集成，不对惩罚后系数使用普通 Rubin 显著性检验。

`model.json` 固定特征顺序、均值、尺度、标准化系数、原尺度系数、中心化常数、目标结局、方向、参考版本、训练指纹与诊断。得分只含影像：`score = X @ raw_coefficients - centering_offset`。临床协变量用于学习时的调整，不混入影像得分。

外部验证检查成员独立、参考与特征顺序，不重新选择权重或惩罚。外部缺失协变量可用已观察验证结局辅助插补，因此输出是固定影像得分的调整后关联，不能冒充未知结局部署预测性能。

主结局复现、方向一致、训练及验证诊断通过后，才标记可进入正式机制注释。最终训练 MI、各训练折 MI 或验证模型存在未解决警报时，结果仍保存但注释资格为否。另保存每套系数与集成系数的方向余弦、逐特征入选频率和符号一致率，不掩盖不稳定权重。

## 6. 哪些内容可改、哪些必须冻结

可在主配置中修改：结局字段与标题、分类等级与缺失编码、协变量集、MI 次数和迭代参数、诊断阈值、纵向开关及优化器、交叉验证折数、alpha、lambda 网格、种子和 FDR 阈值。修改后新建相应统计阶段，不重算不依赖临床配置的病灶特征。

必须保持的分析合同：129 模态的冻结顺序；G/W 含义；观察结局者名单规则；只插补协变量；训练折与验证折隔离；联合 D1/Rubin；失败检验保留；外部验证不调权；临床注释不反向选权。更改这些内容是方法方案变更，应在配置验证层明确拒绝静默切换。

## 7. 方法来源与工程验证

- [mice 插补与 ignore/where 定义](https://amices.org/mice/reference/mice.html)
- [D1 合并检验](https://amices.org/mice/reference/D1.html)；[mitml 参数约束检验](https://search.r-project.org/CRAN/refmans/mitml/html/testConstraints.html)
- [ordinal 模型与说明](https://cran.r-project.org/package=ordinal)
- [ordinalNet 接口](https://search.r-project.org/CRAN/refmans/ordinalNet/html/ordinalNet.html)
- [PSRF 方法：Vehtari 等，2021](https://doi.org/10.1214/20-BA1221)

工程测试入口为 `statistics/selftest.R`、`isolation_test.R`、`longitudinal_test.R`、`fulltest.R` 及 `tests/test_statistics.py`。覆盖独立公式核对 D1/Rubin、四家族分母规则、MI 只改协变量、比例优势诊断、HC3、Jmax、纵向模型、训练隔离、20 组合与 JSON 入口恢复。测试允许缩小 MI、折数与 lambda 数量，但必须将这些工程参数写入测试报告，不能改写正式默认值。
