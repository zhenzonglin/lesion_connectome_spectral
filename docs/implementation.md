# 完整正式项目交付与验收

工程已完成。真实临床分析尚未启动，当前主配置保持待确认状态。所有三例与合成统计结果只用于工程验证，均在 archive；results 当前仅有配置预检。

## 1. 已实现内容

- 冻结 LQT135：135 节点、130 有连接节点、129 非零模态。W 保持原定义；新增 G 和可选 Jmax。
- 完整流程：输入检查 → 病灶特征 → SAS 整合 → 协变量 MI → 逐模态关联 / 20 个组合 → 外部验证 → 空间注释 → 报告。
- 正式代码、config、统计后端、方法说明、健康及注释资源位于主项目；三例与合成验证在 archive。
- 保留精简病例输出、空病灶排除、阶段指纹、监控、安全停止和续跑。没有新增患者病灶质控图或患者人工审核表。

## 2. 实测验收

| 核查 | 实际结果 | 证据路径（项目内） |
|---|---|---|
| Python 回归 | 48 项通过；ID、空病灶、空间、临床整合、指纹、恢复和注释。 | `archive/upgrade/pytest.txt` |
| 三例影像等价 | 既有 W 最大差 0；流线计数逐值一致；G/Jmax 范围与贡献守恒通过。 | `archive/upgrade/imaging_equivalence.json` |
| 输入保护与归档 | 46 个受保护文件校验一致；旧三例与旧检查报告移入 archive。 | `archive/upgrade/migration.json` |
| 真实 SAS 格式接口 | 合成发现/验证各 16 例，实际 sas7bdat 写入、pyreadstat 读取及整合通过。 | `archive/upgrade/sas_roundtrip/verification.json` |
| 统计模块 | 320 例合成资料、258 特征、387 逐模态检验、20 个组合；516/516 个有序拟合的 PO 诊断可计算。 | `archive/statistics_full_final/fulltest.json` |
| 验证与家族 | 48 项预设组合关联；七个月份 mRS、19 项 NIHSS、三个 MoCA 时点及四家族边界通过。 | `archive/statistics_profile/profile_test.json` |
| 防泄漏 | 改变验证侧标签、数值及缺失率，训练插补逐值不变；冻结得分不受验证标签影响。 | `archive/statistics_isolation_final/isolation.json` |
| 纵向与敏感性 | 随机截距、时间交互、保留随机效应的阈值分析和独立 Jmax 分析通过。 | `archive/statistics_longitudinal_5/longitudinal_test.json` |
| 停止与恢复 | 实际正常停止/立即中断均续跑完成；已提交病例不改写，重复启动被拒绝。 | `archive/upgrade/lifecycle.json` |
| 空间注释 | 两个合成空间图各 10000 次 Moran 零模型；286 个计划检验，包括两个缺失 VAChT 位置。 | `archive/annotation_validation/summary.json` |
| 离线路径迁移 | 新路径解包，网络仅 lo；Python/R 自检、注释校验与三例逐值复算全部通过。 | `archive/offline_validation/verification.json` |

工程测试使用缩小的 MI/折数/lambda 网格；正式配置仍为 30 套 MI、5 折、60 个 lambda。测试证明接口和数学/工程约束，不代表全量资料的拟合质量、实际运行时间或临床有效性。

## 3. 资源与存储

健康参考数值残差为 2.37e-15，正交误差为 1.85e-13。全部健康模态已生成 258 张 PNG 与对应脑区表，路径 reference/annotation/healthy。第 73 模态可用作查看示例。

Neurosynth 固定 124 个术语均已从真实元分析图提取。Hansen 预设 19 项，其中 18 项可用；VAChT 缺少受限的 Spreng 原图，保留 NA、原因及检验分母。原始 Schaefer 与 LQT 皮层存在 11343 个体素差异，因此未直接拼接作者分区表。Neurosynth 原分区复现的最小相关为 0.9999999999666。

正式参考目录约 250.00 MiB，含 211.33 MiB 可再生健康图。离线 Python/R 环境压缩包为 748,303,650 字节（713.64 MiB）。三例每例压缩结果与记录约 8.1 KiB，实测每例计算约 4.1—4.5 秒；不外推全量吞吐量。

默认移植包排除病灶、SAS、成员与 ID 映射、患者结果、archive、解包环境和可再生的全部健康图。精确文件清单、包体积及 SHA256 在 portable/lesion_complete_manifest.json；环境包校验与修复说明在 docs/runtime_build.md。

移植实测系统为 Ubuntu 20.04 / x86_64 / glibc 2.31，WSL2 内核 6.6.87.2。测试使用独立网络命名空间，仅有 lo 接口；旧 MATLAB、DSI Studio、FreeSurfer、MRtrix3 路径未加入 PATH。未在其他 Linux 发行版实际运行。

## 4. 用户代码核查入口

1. config.yaml：集中配置与中文分区；先核对路径、ID、SAS 来源和字段。
2. docs/code_review.md：按顺序解释函数、输入、包、处理、输出、可改项。
3. qc/code_checklist.tsv：可编辑的代码审核清单；不对应逐患者病灶审核。
4. docs/clinical.md、docs/statistics.md、docs/annotation.md：临床、模型与空间注释细节。
5. readme.md：启动、监控、停止、续跑和按需导出命令。

## 5. 正式运行前仍需填写

- 实际病灶路径与精确文件模式，或显式 participants.tsv；当前示例 sub-*_lesion.nii.gz 不匹配旧三例特殊文件名。
- 发现与验证 SAS 文件、ID 列/字符串规则及明确的追加或横向合并关系，核查两队列独立。
- 实际变量名、性别/教育/不可测编码，以及各随访 mRS 的 0–6 等级；死亡仅由对应 mRS=6 判定。
- MoCA 字段已经完成教育加分；不再核对版本或重复加分，按配置核对字段是否存在和范围是否符合。

填妥后先运行 `./env/bin/python run.py --run main --preflight`。预检报告会列出缺少的输入，不会把工程验证结果填入临床报告。只有配置和真实资料就位后才启动正式统计。
