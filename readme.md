# LQT135 双通道谐波损伤正式项目

冻结健康参考包含135个解剖节点、130个有连接节点、129个非零模态。每例生成G（局部组织损伤）、W（健康模态连接支撑损失）和可选Jmax。正式模块包括临床整合、协变量多重插补、逐模态关联、20个组合学习、独立验证与空间注释。

三例及合成数据只用于工程验证，保存在archive，不进入发现集或验证集。患者模式形状变化仍为历史探索。不生成患者病灶质控图片或人工审核表。

## 配置与运行

在项目根目录执行。系统为Ubuntu20.04或更新的x86_64 Linux。

```bash
bash setup.sh
./env/bin/python selftest.py --statistics
# 编辑 config.yaml，确认输入路径、ID、SAS映射及临床编码后预检。
./env/bin/python run.py --run main --preflight
./env/bin/python run.py --run main --subjects all --stages preflight,lesion --resume --detach
# 临床配置完成后执行下游。
./env/bin/python run.py --run main --subjects all --stages clinical,imputation,association,patterns,validation,annotation,report --resume --detach
./env/bin/python monitor.py --run main
./env/bin/python stop.py --run main
./env/bin/python run.py --run main --subjects all --resume --detach
```

普通停止在当前病例或统计批次检查点退出；`stop.py --now`只终止本任务核验过的进程。当前双队列续跑显式指定`--subjects all`并保持病例目录不变，各队列子流程会核对原病例清单；新增病例请使用新运行名。修改下游统计参数不会重算病灶；修改上游后，依赖阶段标记过期。旧产物按指纹保留，不混写。

`runtime.workers`控制病灶逐例进程数，范围为1–64；`runtime.threads`控制每个进程的数值库线程数，范围为1–32。当前工作站配置为`workers: 60`、`threads: 1`，用于已分配60个CPU且内存足够的任务；其他机器按可用资源调低。并发改变需要安全停止后使用原运行名续跑，已完成的病灶计算可复用。R统计和组合学习目前仍顺序执行，不会随workers增加并行任务数。

监控默认每30秒持续刷新，显示各阶段及discovery/validation病例计数；任务完成、停止或暂未启动时仍继续观察，按Ctrl+C退出监控。`--once`只显示一次，`--watch`保留兼容，`--interval 30`可覆盖旧配置里的刷新间隔。新任务先发布带运行标识的进度再完成后台启动握手；监控会识别前次运行记录，不把旧失败状态当成本次失败。监控只读，不启动或停止分析。

人数表按发现集、验证集和合计显示“已完成、失败、进行中、待完成、排除、总人数”。已完成仅含成功计算的病例；待完成为尚未开始的排队病例；空病灶计入排除。新运行会保存两队列计划人数，使尚未启动的队列也能显示待完成数量。旧运行或重新核查期间尚无可靠计数时显示“—”，不填零；统计等批量阶段不把模型任务数冒充人数。

续跑预检期间，病灶阶段可能仍指向上轮目录；待开始、过期或受阻阶段不显示该目录中的旧人数。阶段启动后，队列进度早于本次总控启动时间的记录显示“待复核”，直到重新核查后更新；已完成阶段可展示保存的计数。监控显示当前总控PID，`--once`标为“单次查看”。`driver.log`追加保留历次日志，末尾“相同运行正在执行，拒绝重复启动”表示已有总控持锁，应查看当前状态而非重复启动。

`config.yaml`集中管理所有可调参数，中文说明固定定义。代码按序核查见`docs/code_review.md`，SAS合同见`docs/clinical.md`，统计见`docs/statistics.md`，空间对象见`docs/annotation.md`。

正式病灶分为 lesion/discovery 和 lesion/validation 两个目录；两队列可分别配置 pattern、前后缀或命名正则，星号之间完整字符串为ID，保留前导零、大小写及ses文本。验证队列与发现队列ID交集会阻止独立验证。

临床目录按 discovery/validation 分开；程序自动读取各目录全部 .sas7bdat。发现集使用 code_n，验证集使用 c_code+p_code 复合ID；字段映射在 config.yaml 中显式填写。死亡仅由对应 mRS=6 判定，MoCA 已完成教育加分。缺失结局不插补，死后MoCA不填零。

死后仍有 MoCA 分数时默认停止。研究者确认后，可显式设置 `clinical.death.moca_conflict_policy: structural_missing`，仅将首次 mRS=6 当次及之后的 MoCA 设为结构性缺失，保留患者、mRS 和死亡前有效分数。审计表保留处理前的映射分数、首次死亡记录月份及处理动作；死亡后出现非死亡 mRS 的矛盾仍会停止。原始 SAS 不修改。

## 统计与组合

每个结局使用结局已观察、影像完整的固定名单，只对协变量多重插补。默认30套、20次迭代，诊断不足最多100次。同一模态G/W联合D1检验，系数与区间用Rubin规则。有序模型报告比例优势诊断和阈值敏感性；MoCA采用HC3线性模型。纵向mRS为可选有序混合模型，默认关闭。

四个BH家族：NIHSS19项、七次单时点mRS、纵向mRS，以及共用一个家族的3/6/12月MoCA。失败检验保留原因和预设位置，不能通过删失败项缩小分母。

20个组合对应3月mRS及19项NIHSS。全部258个G/W特征进入有序弹性网；临床协变量参与调整且不惩罚。训练折独立标准化及插补，插补不使用结局；默认5折、60个lambda、1SE规则。最终得分只含影像。外部验证冻结权重和标准化参数；验证关联的协变量MI可以用观察结局，不能报告为部署预测性能。

## 结果与脑图

`results/<run>/stages.json`与`report.html`指向`stages/<stage>/<fingerprint>/`产物。每例仅压缩NPZ和JSON记录，包含1375条固定健康边中断计数、g/w/jmax、135脑区覆盖率。汇总列为`g_001…g_129`、`w_001…w_129`。`summary.tsv`中的原始体积单位为mm³。

```bash
./env/bin/python export.py --healthy --mode 73
./env/bin/python export.py --healthy --squared --mode 73 --name healthy_square_73
./env/bin/python export.py --run main --subject 0001 --mode 73 --channel g
./env/bin/python export.py --run main --subject 0001 --mode 73 --channel w
```

健康图是有符号载荷或平方载荷。患者及组合图是损伤贡献，不能用健康特征向量线性相加替代。组合图展示外部队列对应结局已观察患者的平均G/W贡献；两图节点和等于平均影像得分加冻结中心化常数。零度节点为NaN。首版注释推断仅在100个皮层节点进行；小脑、脑干和皮层下仍参加损伤计算。

真实注释资源的来源、许可、校验和可用性见`reference/annotation/manifest.json`。受限资源保留NA与计划检验数，不以别的分区均值冒充精确重提取。空间注释不反向选择临床权重。

## 离线移植

```bash
./env/bin/python package.py --name lesion_complete
```

移植包包含代码、Python/R环境、精简健康索引、注释资源和说明，默认不含病灶、SAS、ID映射、患者结果或archive。解包到最终路径，执行setup和自检，原始数据另行转移。环境重建参考environment.yml，build_runtime.py制作conda-pack环境。正式计算无需MATLAB、DSI Studio、FSL、FreeSurfer或MRtrix3。

本轮实测验收见`docs/implementation.md`，详细测试证据在archive；临床发现必须等待真实资料配置后另起运行。
