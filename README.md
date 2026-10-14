# 工作站轻量更新分支

本分支使用独立历史，只包含 5 个代码更新文件、SHA256 清单、维护脚本和测试；没有主分支历史、参考图谱、环境、病灶、SAS 或患者结果。适用于已经部署的 lesion_connectome_spectral 项目，不能独立运行分析。

代码来源：`7cf4af23aa030a114235a1d847a26956f82acb9e`。内容包括病灶调度修复、60 工作进程支持、30 秒持续监控、分队列人数表，续跑时不把上轮计数显示为本轮状态，以及显式启用并审计死后 MoCA 的结构性缺失处理。

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

安装前核查停止状态、补丁 SHA256 和 Python 语法；旧代码备份到 `archive/workstation_updates/`，再替换五个白名单文件。异常时恢复已替换文件。配置、病灶、SAS、参考资源、环境和已有结果均保留。`functions/cohort.py` 不在补丁中，工作站已做的递归扫描修改会保留。

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

该选项将首次 mRS=6 当次及之后的非缺失 MoCA 设为结构性缺失，保留患者、全部 mRS 和死亡前有效分数，不修改原始 SAS。`audit/death_conflicts.tsv` 记录处理前的映射分数、首次记录死亡的月份、处理动作和策略；日志输出涉及人数及观测条数。死亡后出现非死亡 mRS 的矛盾仍停止。其他可用结局不因这项处理而删除，结局不插补。

安装代码并设置策略后，可使用原运行名仅续跑临床及下游，病灶输入和代码一致时复用病灶结果：

```bash
./env/bin/python run.py --config "$LCS_CONFIG" --run "$LCS_RUN" \
  --subjects all \
  --stages clinical,imputation,association,patterns,validation,annotation,report \
  --resume --detach
```

## 开发验证

在现有项目根目录执行 `./env/bin/python -m pytest "$LCS_PATCH/tests" -q`。测试覆盖进程范围、子进程退出、残留 worker、运行锁、安装备份、SHA256、失败回滚及图像头检查。测试使用临时项目和模拟进程，不能替代工作站的真实数据检查。
