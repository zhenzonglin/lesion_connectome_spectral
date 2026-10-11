# 空间对象与注释模块

正式入口为 `functions/annotation.py`。它接收已冻结的健康参考和临床组合定义，不参与临床组合的选择。`annotate.py` 提供按需脑图导出。原始注释图及一次性元分析构建环境位于 `archive/annotation_source`，移植和正式运行只需要 `reference/annotation` 中的紧凑脑区表、来源与校验值。

## 1. 三类空间对象

|对象|函数|输入|输出与解释|
|---|---|---|---|
|健康背景|`export_healthy_modes`|冻结 `reference/modes.npz` 中的129个非零模态|每个模态的有符号载荷和平方载荷，各有 `regions.tsv` 与 `map.png`；正式资源位于 `reference/annotation/healthy`|
|患者节点损伤贡献|`loss_contributions`|135维覆盖率、1375条固定边的中断计数|135×129 的 G 与 W 贡献，仅按需计算；排除节点为NaN|
|经验证组合贡献|`combination_maps`、`mean_combination_maps`|患者贡献、冻结258维系数及训练均值/标准差|G和W各一张135区贡献图；按指定患者集合逐例平均，不储存巨大三维贡献数组|

G节点贡献为 `rho_i × u_ik²`。W节点贡献为每条受损边的健康模态能量损失各分一半给两个端点。G/W节点值分别求和，必须恢复该患者该模态的G/W损伤量。健康模态整体换符号不会改变任何损伤贡献。

设标准化G/W特征合并为258维，冻结系数为 `a`，训练均值为 `mu`，标准差为 `sigma`。原尺度系数为 `a/sigma`。各模态的节点损伤贡献乘原尺度系数后再求和，即得到组合贡献图。

两张图的节点总和减去固定的 `sum(a × mu / sigma)`，必须等于患者影像组合得分。临床协变量不进入这张影像得分图。不能用 `sum(a_k × u_k)` 代替损伤贡献图。

队列平均贡献图同时受冻结组合系数和该队列的病灶分布影响。它不是组合系数唯一决定的健康脑图，也不是病灶造成局部神经活动升降的测量。G/W分图展示正负贡献；5个未进入健康分解的节点保持缺失。

## 2. 分区一致性与原始图重提取

Hansen仓库发布的100区表使用 Schaefer100 / 7Networks。项目100个皮层标签的名称、左右顺序和CBIG定义一致，但项目为了让深部灰质和皮层互斥，11343个原Schaefer体素被深部标签覆盖。与原分区的Dice均值约0.9943、最低约0.8764。因此，不能直接把已发表的100行数据当作本项目皮层体素的均值。

`prepare_annotations.py::alignment_audit` 检查标签名称、脑区顺序、仿射和实际体素。`parcel_means` 将项目标签以最近邻重采样到原始注释图的网格，然后计算各区均值，不插值注释信号。PET只使用有限且非零的覆盖体素；Neurosynth的零z值是有效数据，保留。末尾单元素维度的PET图先压缩为三维。

Hansen原图按官方 [`make_receptor_matrix.py`](https://github.com/netneurolab/hansen_receptors/blob/main/code/make_receptor_matrix.py) 聚合为固定19项。多来源图先z标准化再按原源码权重组合。必须保留原始PET研究引用，项目保存来源表和Hansen许可文件（CC BY-NC-SA 4.0）。

作者VAChT源码实际为 `(3*z3+4*z4+z5+z18)/30`，与“按样本数加权”的文字不完全一致。本项目主表忠实复现实际源码；另存 `(3*z3+4*z4+5*z5+18*z18)/30` 敏感性表，不加入主19项。Spreng的3人VAChT原始体积图当前受上游OSF权限限制；若不能取得，整项VAChT保留NA及`resource_unavailable`，不拿原分区均值替代，也不把剩余三张原图改组成新的主VAChT指标。该项仍保留在预设多重比较的分母中。

Neurosynth术语集直接冻结Hansen发布CSV的列名。README声称123项，但实际CSV是124项；不任意删一列凑数。利用作者课程提供的[原始流程](https://github.com/netneurolab/ipn-summer-school/blob/main/lectures/2021-07-02/13-15/fetch_and_parcellate_neurosynth.py)，调用Neurosynth 0.3.8、固定2018数据库 `e8f27c4`、6mm球、术语频率阈值0.001，生成未阈值化 `association-test_z`，再按项目真实皮层分区提取。现代NumPy/Pandas/NiBabel的兼容别名仅处理弃用API；不改变统计公式。数据库、原图及紧凑表均有SHA256。

## 3. 空间零模型和两个注释家族

主推断只使用全部100个匹配皮层节点。皮层下、小脑、脑干保留覆盖说明，不混入皮层统计。任一皮层节点缺失时阻止该图推断，不静默删脑区。

默认 `hemisphere_singleton_moran` 使用半球内脑区MNI质心欧氏距离的倒数作为对称空间权重。通过Moran特征向量的随机符号翻转生成10000张空间零图。算法遵循[BrainSpace实现](https://github.com/MICA-MNI/BrainSpace/blob/master/brainspace/null_models/moran.py)与Wagner和Dray的singleton定义，使用SciPy双精度无常数子空间，避免额外VTK环境。

每张零图保留每侧半球均值、方差和该权重定义下的Moran空间二次型。程序逐图核查保真误差并保存诊断。它不是球面旋转，也不代表皮层测地距离；不宣称保留原图的全部经验分布。该距离选择与误差会写入每次报告。增加另一种零模型需要新方法版本，不能仅修改名字而继续执行旧算法。

计算Pearson空间相关，双侧空间p值为 `(1 + 零图中绝对相关不小于观测的数量)/(1+零图数)`。不得使用把100个脑区当独立样本的普通相关p值。

Neurosynth和Hansen分别做BH-FDR，覆盖各自预设的所有组合、G/W通道和术语/受体。资源不可用、常量脑图和未复现模式均记录NA及原因，仍保留预设检验位置，避免通过删除失败项减小校正分母。

## 4. 阶段输入输出合同

`run_annotation(root, config, input_dir, output_dir)` 读取 `input_dir/maps_manifest.json`。其中 `maps` 的每项必须包含 `pattern_id`、`channel`（g或w）、`map_file` 和 `validated`。TSV必须有完整唯一的1至135 `index` 和 `value`。只能引用该输入目录内文件；节点顺序依编号连接。

仅`validated=true`的组合进入推断。上游应依据外部验证的方向、家族校正和诊断产生此状态，不能因一张图好看而手工升级。输出包括 `associations.tsv`、展示PNG、`null_diagnostics.json`及 `summary.json`。CSV/TSV的完整结果优先于只展示绝对相关较大项目的概览图。

以下参数来自 `config.yaml` 的 `annotation` 分区：资源目录、零模型名称、随机次数、随机种子、图像DPI和色图。资源的脑区数量、术语名单、原数据版本与SHA属于冻结资源校验属性，不能随意修改这些数字改变分析范围。

```bash
# 核验冻结注释资源，不联网
./env/bin/python prepare_annotations.py --verify

# 按需展示患者第73模态的G/W贡献，输出新目录
./env/bin/python annotate.py --patient results/myrun/subjects/0001/result.npz \
  --mode 73 --output results/exports/patient_0001_mode73

# 已验证组合的正式注释，通常由run.py自动调用
./env/bin/python annotate.py --input results/myrun/combination_maps \
  --output results/myrun/annotation
```

首次开发的三例与合成测试位于 `archive`，不进入任何正式发现或验证队列。患者贡献导出不是病灶质控，也不会批量生成患者病灶叠加图。
