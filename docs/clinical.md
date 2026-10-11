# 2026 双队列正式输入合同

正式配置自动扫描 clinical/discovery 和 clinical/validation 中的全部 .sas7bdat，不需要逐个填写 SAS 文件名。每个队列在 clinical.cohorts 中定义身份列：发现集使用 id.field: code_n；验证集使用 id.fields: [c_code, p_code] 与 separator: "" 拼接复合 ID。字段名和大小写必须与真实 SAS 列一致，代码不会根据相似名称自动别名。

病灶路径由 paths.lesion.discovery、paths.lesion.validation 指定，文件名规则由 subjects.cohorts.<cohort>.pattern 与 id_rule 指定。每个队列先独立扫描、去重并和同队列临床 ID 匹配；发现/验证 ID 交集会阻止验证。空文件和全零掩膜记录为 excluded_empty_lesion，不生成零损伤特征。

验证队列字段覆盖写在 clinical.cohorts.validation.field_map。目前配置填写了共同基础变量和 m3_mrs36/m6_mrs36/m12_mrs36；缺少的 NIHSS 分项、MoCA 或长期 mRS 会进入共同变量审计表，不能静默映射到发现集字段。只有两队列均存在且编码合同一致的结局进入独立验证。死亡仅由对应访视 mRS 等于 6 判定；不读取独立 death 字段。MoCA 已完成教育加分，不再检查版本或重复加分。

旧运行的来源操作说明已迁入 archive，不作为当前正式字段合同。当前正式配置只使用本页明确的双队列字段和 mRS=6-only 死亡规则。

# 临床输入与字段核查

本模块只导入、核对、整理临床记录。原始 SAS、病灶和历史结果只读；临床数据不得进入默认离线包。主配置默认没有真实 SAS 文件，两个队列均为 `confirmed: false`。临床阶段会明确列出缺项；这不是执行失败后自动改用测试数据。

## 1. 需要填写的配置

正式流程会自动扫描两个目录中的所有 SAS 文件：

- clinical/discovery/*.sas7bdat：发现集，ID 默认取 code_n。
- clinical/validation/*.sas7bdat：验证集，ID 默认按 c_code 与 p_code 无分隔拼接。

因此不需要在配置中逐个填写 SAS 文件名。只需把文件放入对应目录，并在 config.yaml 的 clinical.cohorts.<cohort>.id 中核对 ID 列。若实际 ID 不是默认形式，明确填写 field 或 fields、separator、kind 和必要的 width。代码不会根据相似字段名自动别名。

验证集只使用两队列都存在且编码合同一致的变量。当前配置把验证集的 3、6、12 月 mRS 映射到 m3_mrs36、m6_mrs36、m12_mrs36；缺少的 NIHSS 分项、MoCA 和长期 mRS 会记录在共同变量审计中，不生成伪造的缺失结局或零值。

如果同一队列有多张 SAS 表，程序按 subject_id 横向合并；同名字段出现非缺失冲突会停止并写入审计表。每张表内重复 ID、发现/验证队列交集和临床—影像不匹配都会明确记录。成员清单是可选的 TSV，只能有唯一 subject_id 列。

## 2. 变量与死亡规则

`clinical.variables` 将内部名称、标题、SAS 字段、类型、合法范围/等级、缺失编码分开。候选映射来自用户变量表，不代表实际 SAS 已核实。`continuous` 的 `range` 包含两端；`binary/nominal/ordinal` 必须明确 `levels`。支持字符串分类等级。外部队列只允许覆盖 `field/source/missing/recode`，不能改变同一个内部变量的含义或类型。

`pyreadstat.read_sas7bdat(user_missing=True)` 先保留 SAS 特殊缺失信息，再根据元信息转为空值；元信息、标签与来源 SHA256 写入清单。`9` 仅在 NIHSS 5A、5B、6A、6B、10 中作为不可测编码；其他 NIHSS 字段出现 `9` 会报错，不静默转缺失。EDUC 的 `98` 为缺失；教育类别不换算成年数。

`clinical.derived` 仅允许：

- `log1p`：对原始病灶体积 `lesion_volume_mm3` 计算自然对数 `log(1+体积)`，单位仍明确为 mm³ 输入。
- `sum_complete`：组成项都有效才求和；任一不可测或缺失则合并项缺失。

入院 NIHSS 15 分项和四个合并项是独立结局。mRS 调整入院 NIHSS 总分和卒中前 mRS；MoCA 只调整年龄、性别、对数体积和教育。

死亡状态只由对应随访的 mRS 是否等于 6 推断，不读取独立 death 字段或另行配置死亡编码。记录 death_<month>m 状态；mRS=6 后再次出现非缺失的非 6 mRS 会被报告为冲突。mRS 缺失不自动补成 6，也不因缺失推断死亡。

死亡后缺失 MoCA 额外标记 `moca_*m_structural_missing`，不插补、不记零。死后非缺失 MoCA 报错。MoCA 已在上游完成教育加分；本项目不识别版本、不重复加分。死亡规则只对启用分析所需月份执行。

## 3. 逐步代码核查表

| 顺序、函数 | 输入 | 包与处理 | 输出 | 固定内容 / 可修改内容 |
|---|---|---|---|---|
| `configuration.load_config` | YAML 文件 | PyYAML 安全读取；检查分区、双队列规则、类型、四家族、20 目标 | Python 配置字典 | 不执行表达式；变量映射和参数可改；健康数量为冻结属性 |
| clinical.preflight_clinical | 项目根、配置、队列名 | 标准库核查 SAS 路径、ID、mRS=6-only 死亡规则及死后 MoCA 冲突 | 缺项字符串列表 | 实际路径、字段和量表编码需人工核对 |
| `clinical.normalize_ids` | SAS ID 列、ID 规则 | pandas；显式整数补零或对应表；保持文本身份 | `subject_id` | 不猜测、不去 ses；补零宽度和映射可配 |
| `clinical.read_sas` | SAS 路径 | pyreadstat；保留标签和特殊缺失元信息 | 数据框、来源元信息 | 不使用显示标签覆盖实际编码 |
| `clinical.load_sources` | 配置的 SAS 列表 | pandas；base、append、join；来源命名空间；重复核查 | 合并原始宽表、来源清单、未匹配/重复表 | 重复不取第一条；来源操作必须配置 |
| `clinical.map_variables` | 合并宽表、变量定义 | pandas/NumPy；字段映射、缺失编码、可选 recode、类型与范围核查 | 内部名称宽表、编码异常表 | 异常不能静默处理；范围/等级按量表核实 |
| clinical.check_death | 映射表、各随访 mRS 与 MoCA | pandas；mRS=6-only 死亡状态、跨月份一致性和死后 MoCA 冲突核查 | 死亡状态、结构性缺失标记、冲突表 | 不填充结局；mRS 缺失不推断死亡 |
| `clinical.derive_variables` | 内部变量、原始病灶体积 | NumPy `log1p` 与 pandas 完整项求和 | 四个 NIHSS 合并项、对数体积 | 支持运算固定；派生输入在 config 指定 |
| `clinical.prepare_clinical` | 以上配置、影像特征目录、输出目录 | pandas 一对一连接；特征完整与范围核查；各结局观察数统计 | `analysis.tsv.gz`、`manifest.json`、`audit/*.tsv` | 不插补结局/影像/体积；不生成病灶质控图片 |
| `clinical.assert_disjoint` | 发现/验证原始成员集合 | Python 集合求交 | 无重叠通过，否则报错 | 不能关闭后假称独立验证 |

函数均位于 `functions/clinical.py`；配置读取位于 `functions/configuration.py`。依赖和方法依据：[pyreadstat 官方代码](https://github.com/Roche/pyreadstat)。统计推断由后续 R 模块实现，本模块没有使用合成病例估计真实效应。

## 4. Python / R 接口和输出位置

```python
config = load_config(root / "config.yaml")
issues = preflight_clinical(root, config, cohort="discovery")
manifest = prepare_clinical(root, config, feature_dir, output_dir, cohort="discovery")
```

`feature_dir` 中必须有 `mode_features.tsv.gz` 和 `summary.tsv`。前者包含 `subject_id,g_001...g_129,w_001...w_129`；后者包含 `subject_id,lesion_volume_mm3`。缺失、非有限或范围外的影像特征会报错。全零病灶已经在影像阶段排除，不允许被当成正常的零特征。

输出全部写在调用者传入的 `output_dir`，正式总控放在本次 `results/<run>/...`；测试传入 `archive/...`。`analysis.tsv.gz` 只保留影像成功且临床匹配患者，包含原始缺失值；后续每个结局独立筛选观察者。`manifest.json` 包含全部原始队列成员、影像匹配成员、来源 SHA256、258 个固定特征顺序和 R 所需结局定义。R 接口的变量 `type/levels` 由字典提供，类别不会被误当连续量。

核查表包括来源未匹配、影像临床未匹配、编码异常、死亡冲突、缺失分布、逐结局样本流向。若存在重复或字段冲突，先写相关表再明确停止。导入阶段不会插补缺失协变量；全部缺失的必要协变量会由后续统计预检阻止该结局拟合。

默认四个家族是 NIHSS 19项、mRS七时点、纵向mRS 和 MoCA 三个时点；三个 MoCA 时点共用一个 FDR 家族，纵向默认关闭。`outcome_definitions(config)` 将家族和协变量集展开成 R 所需列表。`clinical_input_files(root, config)` 返回所有 SAS、成员表及 ID 对应表，供总控计算阶段指纹。
