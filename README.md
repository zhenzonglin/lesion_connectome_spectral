# 工作站轻量更新分支

本分支使用独立历史，只包含 6 个代码更新文件、SHA256 清单、维护脚本和测试；没有主分支历史、参考图谱、环境、病灶、SAS 或患者结果。适用于已经部署的 lesion_connectome_spectral 项目，不能独立运行分析。

代码来源：`022458814a3d673caa7ba82769bc8384a19bb94b`。内容包括病灶调度修复、60 工作进程支持、30 秒持续监控、分队列人数表，续跑时不把上轮计数显示为本轮状态，以及显式启用并审计死后 MoCA 的结构性缺失处理和按队列排除 mRS 矛盾患者的 mRS/MoCA 结局。

## 仅导出汇总审阅包（2026-10-18）

若分析已经完成、只需要将关键结果提供给审阅者，使用独立 `export_review.py`。
它从当前完成阶段中重建白名单汇总表和健康参考，不打包原始临床资料、患者ID、逐人结果、RDS或原HTML。

```bash
cd /data/usersdir/linzhenzong/lesion_connectome_spectral || exit 1
LCS_PATCH=/data/usersdir/linzhenzong/lesion_connectome_spectral_hotfix_20261014
git -C "$LCS_PATCH" pull --ff-only origin hotfix/workstation-lite || exit 1
./env/bin/python "$LCS_PATCH/export_review.py" \
  --project "$PWD" --config config.workstation.yaml --run main_20261014 \
  --output "$PWD/review_export_$(date +%Y%m%d_%H%M%S).zip"
```

加 `--dry-run` 可先预览。已有同名ZIP不会覆盖。此脚本直接运行，无需执行下方的 stop/install/check/run，也不向完整工作站项目复制补丁。
字段合同、隐私排除项和缺失说明见 [审阅导出文档](docs/review_export.md)。
导出器与测试为独立新增文件；下方原6个payload及其清单保持不变。

## 一次完整的停止、检查与续跑

以下命令在同一个 Bash 终端依次执行，每步成功后再继续。本次保留原运行名和所有已有结果；更新和检查期间不要从另一个终端启动分析。

### 1. 下载轻量分支

```bash
export LCS_ROOT=/data/usersdir/linzhenzong/lesion_connectome_spectral
export LCS_PATCH=/data/usersdir/linzhenzong/lesion_connectome_spectral_hotfix_20261014
export LCS_CONFIG=config.workstation.yaml
export LCS_RUN=main_20261014

git clone --depth 1 --single-branch --branch hotfix/workstation-lite \
  https://github.com/zhenzonglin/lesion_connectome_spectral.git "$LCS_PATCH"

cd "$LCS_ROOT"
```

若补丁目录已经是本分支的干净克隆，后续更新使用 `git -C "$LCS_PATCH" pull --ff-only`。不要在完整项目目录切换到这个精简分支。

### 2. 停止本用户在本项目中的全部进程

```bash
./env/bin/python "$LCS_PATCH/workstation.py" stop \
  --project "$LCS_ROOT" --config "$LCS_CONFIG"
```

脚本识别本项目的总控、队列、后端、监控和工作进程，包括当前运行名以外的同项目任务。先向协调进程发送终止请求，等待 30 秒；仍存活的已核验进程会被终止。每次发送信号都复核 PID 创建时间及当前用户，不能只按 `python`、`R` 或路径子串匹配。

成功必须显示“活动进程 0，持有的运行锁 0”。若仍有进程或锁，命令返回非零，不要继续。不会删除运行锁文件或已有结果。未提交病例和被中断的统计阶段将在续跑时重新计算。

### 3. 安装代码补丁

```bash
./env/bin/python "$LCS_PATCH/workstation.py" install \
  --project "$LCS_ROOT" --config "$LCS_CONFIG"
```

安装前核查停止状态、补丁 SHA256、Python 和 R 语法；旧代码备份到 `archive/workstation_updates/`，再替换六个白名单文件，包括 `functions/clinical.py` 和 `statistics/imputation.R`。异常时恢复已替换文件。配置、病灶、SAS、参考资源、环境和已有结果均保留。`functions/cohort.py` 不在补丁中，工作站已做的递归扫描修改会保留。

### 4. 启动前检查

```bash
./env/bin/python "$LCS_PATCH/workstation.py" check \
  --project "$LCS_ROOT" --config "$LCS_CONFIG" --workers 60
```

检查当前进程与运行锁、安装代码、CPU 亲和范围、`runtime.threads: 1`、冻结参考及数值定义、两队列文件清单与 ID、SAS 必需字段、R 依赖，以及全部病灶图像头中的三维形状、有效仿射矩阵和 1 mm 体素。扫描文件并计算校验值可能需要时间，不启动病灶计算。

只有出现 `CHECK PASSED` 且命令正常退出后再启动。检查记录保存在 `archive/workstation_checks/`。零字节病灶按项目原规则记录为后续排除对象。图像头检查不验证解剖配准，也不读取全部体素检查二值性；正式病灶阶段保留原有检查。内存需求依病例和工作进程而变，检查输出可用内存但不把它当作内存足够的保证。

### 5. 用原运行名续跑所有已启用阶段

```bash
./env/bin/python run.py \
  --config "$LCS_CONFIG" --run "$LCS_RUN" \
  --subjects all --workers 60 --stages all --resume --detach

./env/bin/python monitor.py \
  --config "$LCS_CONFIG" --run "$LCS_RUN" --interval 30
```

启动只执行一次；“相同运行正在执行”表示已有总控持锁。`--stages all`遵循配置中的阶段列表和各模块启用开关。仅指纹一致且完整的结果会复用；修正病灶文件或路径可能生成新的病灶阶段目录，包括重新计算发现集，旧结果保留。60 工作进程用于病灶逐例任务，R 后端仍按原方式执行。

监控每 30 秒刷新，显示已完成、失败、进行中、待完成、排除和总人数；等待本轮核查时显示 `—`。按 Ctrl+C 退出监控，后台分析继续。

## 单独查看是否已经全部停止

```bash
./env/bin/python "$LCS_PATCH/workstation.py" status \
  --project "$LCS_ROOT" --config "$LCS_CONFIG"
```

## 已确认的死后 MoCA 处理

临床更新支持 `clinical.death.moca_conflict_policy`。默认或 `error` 保持严格拦截。研究者明确确认按 mRS=6 处理后，在现有配置的 `clinical.death` 中加入：

```yaml
moca_conflict_policy: structural_missing
```

该选项将首次 mRS=6 当次及之后的非缺失 MoCA 设为结构性缺失，保留患者、全部 mRS 和死亡前有效分数，不修改原始 SAS。`audit/death_conflicts.tsv` 记录处理前的映射分数、首次记录死亡的月份、处理动作和策略；日志输出涉及人数及观测条数。死亡后出现非死亡 mRS 的矛盾默认仍停止；下面的结局排除是独立选项。其他可用结局不因 MoCA 缺失处理而删除，结局不插补。

安装代码并设置策略后，可使用原运行名仅续跑临床及下游，病灶输入和代码一致时复用病灶结果：

```bash
./env/bin/python run.py --config "$LCS_CONFIG" --run "$LCS_RUN" \
  --subjects all \
  --stages clinical,imputation,association,patterns,validation,annotation,report \
  --resume --detach
```

## 验证集 mRS 矛盾患者的结局排除

研究者选择：存在跨月 mRS 死亡矛盾的验证集患者退出全部 mRS 和 MoCA 结局分析，保留其他分析及审计。不要删除患者的临床行或病灶文件，也不要将后续 mRS 回填为 6。默认仍为严格拦截；在工作站现有 `config.workstation.yaml` 的 `clinical.cohorts.validation` 下加入以下 `death` 子项：

```yaml
clinical:
  cohorts:
    validation:
      death:
        mrs_conflict_policy: exclude_mrs_moca
```

这是定位示例，请合并到已有分区，不新增重复的 `clinical`/`validation`。保留验证集原来的 ID 和字段映射，也保留已确认的 `clinical.death.moca_conflict_policy: structural_missing`。发现集默认仍拦截 mRS 矛盾。

更新顺序：先 `git -C "$LCS_PATCH" pull --ff-only origin hotfix/workstation-lite`，再执行本页的 `stop`、`install`，修改配置并续跑 clinical 及下游。安装同时更新 Python 标记和 R 筛选；只设置 YAML 而未更新这两个文件不会实现排除。

患者及处理前分数仍保存在分析宽表；`mrs_moca_outcome_excluded=1` 在 R 的结局入口排除所有 mRS/MoCA 访视，包括早期访视和纵向 mRS。筛选发生在插补、拟合、组合学习和外部验证之前。NIHSS 等其他结局不使用这个排除标记。矛盾患者的推断死亡标记置空，不将其当成已核实死亡。

完成 clinical 后，在项目根目录执行下列只读核验；始终使用本轮阶段指针，不从旧目录取审计：

```bash
./env/bin/python - <<'PY'
import json
from pathlib import Path
import pandas as pd
root = Path.cwd()
state = json.loads((root / 'results/main_20261014/stages.json').read_text())
record = state['stages']['clinical']
print('本轮 clinical 状态:', record['status'])
if record['status'] != 'completed':
    raise SystemExit('等待 clinical 完成；若失败，请查看本轮 execution.log。')
folder = root / record['path'] / 'validation'
exclusions = pd.read_csv(folder / 'audit/outcome_exclusions.tsv', sep='\t', dtype={'subject_id': str})
print('临床来源中排除人数:', exclusions.subject_id.nunique())
print('涉及 mRS 冲突条数:', int(exclusions.conflict_records.sum()))
print('排除清单:', folder / 'audit/outcome_exclusions.tsv')
flow = pd.read_csv(folder / 'audit/sample_flow.tsv', sep='\t')
print(flow[['outcome', 'image_matched', 'excluded_outcome', 'observed_outcome', 'missing_outcome']].to_string(index=False))
PY
```

针对已报告的这批记录，临床来源预期为 23 人、25 条冲突；实际以本轮审计为准。`excluded_outcome` 只统计有有效影像并匹配到临床的患者，可能少于来源排除人数。排除规则不硬编码 23，也不据此推定最终模型样本量。发现集先前 7 人的 8 条死后 MoCA 结构性缺失策略保持有效。

## 开发验证

在现有项目根目录执行 `./env/bin/python -m pytest "$LCS_PATCH/tests" -q`。测试覆盖进程范围、子进程退出、残留 worker、运行锁、安装备份、SHA256、R 语法拒绝、失败回滚及图像头检查。测试使用临时项目和模拟进程，不能替代工作站的真实数据检查。
