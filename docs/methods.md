# 方法、符号与来源

## 冻结参考

LQT 原版仓库提交 `050bc257d686d7e7cba28262429e583bf5be7454`，使用 HCP-842 群体平均扩散流线图谱。它不是 842 个独立被试 SC 的平均。LQT135 由 100 个皮层、8 个深部灰质、26 个小脑和 1 个整体脑干分区构成。深部标签已根据实际图像核对，原标签仍保留在 `original_label`。

健康图 135 节点、1375 条无向边。130 个有连接节点参与健康归一化拉普拉斯分解，得到 129 个非零模态。5 个零度脑区为 Cerebelum_3_L、Cerebelum_10_L、Vermis_1_2、Vermis_9、Vermis_10；不补连边，也不将其缺失特征报告为零。

端点规则是 `end`；经过脑区但未在其中终止不算该脑区端点。主图不增加皮层局部邻接，不二值化，不做连通阈值筛选。图谱外病灶和未匹配端点独立报告。

## 连接和谱

患者编号 p，脑区编号 i、j，模态编号 k。健康流线计数 n_ij；中断计数 m_pij；固定脑区体积 v_i（mm³）。

`c0_ij = n_ij / (v_i + v_j)`；`cp_ij = (n_ij - m_pij) / (v_i + v_j)`。

健康连接强度 d0_i 是 c0 每行之和，d0 是其对角矩阵。健康拉普拉斯 `l0 = I - d0^(-1/2) c0 d0^(-1/2)`，`l0 u_k = lambda_k u_k`，每个 u_k 为单位长度向量。

`delta_c = c0 - cp`；`delta_k = d0^(-1/2) [diag(delta_c 1) - delta_c] d0^(-1/2)`。

逐模态支撑损失 `s_pk = u_k^T delta_k u_k / lambda_k`。固定健康强度、健康模态和特征值；不重新归一化患者图。它等于被删除连接在该健康模态图能量中所占比例，是本项目扩展，不是 Preti 的 fMRI SDI，也不是神经活动或功能丧失比例。

每条无向边能量 `delta_c_ij (u_ki/sqrt(d0_i) - u_kj/sqrt(d0_j))²`，除以 lambda_k 后均分给两个端点。得到网络贡献 w_pik；对所有有效脑区求和等于 W_pk（原变量 support_loss / s_pk）。按需导出的 NIfTI 使用每脑区常值，空间注释直接采用脑区值，避免体积重复加权。

新增局部通道 `G_pk=sum_i rho_pi*u_ik²`；节点贡献 `g_pik=rho_pi*u_ik²`。G和W均在0—1，但不是可以互换的损伤机制。`Jmax_pk=max(G_pk,W_pk)`是可选描述/敏感性特征，主模型只同时纳入同一模态G和W。三个量都是操作性模型定义，不是直接测得的功能缺损。

组合得分 `T_p=sum_k a_k*(G_pk-mu_Gk)/sigma_Gk + b_k*(W_pk-mu_Wk)/sigma_Wk`。a、b为标准化特征权重。空间图按`a_k/sigma_Gk`及`b_k/sigma_Wk`分别加权节点G/W贡献，两图节点总和为T加固定中心化常数。组合权重允许为负，空间图表示有符号统计贡献。临床展示对象为外部队列对应结局已观察患者的平均贡献，图注报告队列和样本数。

灰质覆盖率 rho_pi 是病灶与脑区交集体积除以脑区体积，独立保留，不乘入 SC。整体脑干包含灰白质混合组织，其覆盖率单列。路径命中使用完整病灶，不限定病灶组织类别。

全谱汇总 `sum(u_k^T delta_k u_k)/sum(lambda_k)` 等于 130 个分析节点的平均相对连接强度损失；逐模态曲线才反映不同健康空间模式的差异。

## 坐标和精简索引

原生 DSI 点坐标是 TRK 原始坐标除以 voxel_size，再乘 FIB 仿射；FIB 转置 trans 与 QA NIfTI 仿射相符。端点在原生 QA 网格使用 floor(x+0.5) 分配。索引通过 `inv(template.affine) @ qa.affine` 转换到固定 MNI 模板网格。

完整线段 DDA 与既有 `geometry.segment_hit` 采用同一半体素边界、1e-12 交点容差及同时过面规则。先去重体素—流线组合键，再由全部病灶体素的命中流线取并集；同一条流线反复穿越只计一次。低于两个采样点的 37 条流线不参与路径或端点计算。

138445 条流线：37 条过短、38686 条至少一端未匹配、8170 条同脑区、91552 条有效脑区间流线。端点原始计数已与 DSI 原版输出逐元素核对。索引加速属于本项目实现，不能声称运行时仍调用了原版 LQT。

来源文件中的旧绝对路径仅作为历史来源记录，不用于运行。模板版本声明与物理网格核查保留；图像仍需人工检查。三例幕上病灶和模拟幕下对照不构成真实小脑/脑干患者临床验证。历史 Gibbs 模板版本限制保留在旧归档，不转用为当前 HCP-842 的元信息结论。

## 来源

- LQT：https://doi.org/10.1016/j.nicl.2021.102639
- LQT 代码：https://github.com/jcgriffis/lesion_quantification_toolkit/tree/050bc257d686d7e7cba28262429e583bf5be7454
- HCP-842：https://doi.org/10.1016/j.neuroimage.2018.05.027
- Preti 健康谐波：https://doi.org/10.1038/s41467-019-12765-7
- Preti 作者代码：https://github.com/gpreti/GSP_StructuralDecouplingIndex
- NiBabel 重采样：https://nipy.org/nibabel/reference/nibabel.processing.html
- 离线环境：https://conda.github.io/conda-pack/

LQT README 的 GPL-3.0-or-later 声明和 GPL 原文随 reference 保存。原始来源的图谱与数据条款仍适用；本项目不改变其授权。Python、Conda 依赖的许可证随环境及其 conda-meta 记录保留。

## 用户确认的实施修正

空病灶状态为 excluded_empty_lesion，保留病例记录并从分析特征表排除，不生成零值 SC 或零值扰动。取消患者病灶质控图及人工审核，仅保留数学计算必需的格式与坐标检查。既有质控图和人工意见在历史归档中保留。
