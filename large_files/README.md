# 大型离线资源（不进入普通 Git）

本目录保存可选的大型二进制快照。它们不包含真实患者病灶或 SAS 文件，但单个文件超过 GitHub 普通 Git 的 100 MiB 限制，因此不提交到仓库。每个文件的用途、大小和 SHA256 见 manifest.tsv。

- portable/lesion_complete.tar.gz：不含患者数据的完整移植快照，包含代码、冻结健康参考、文档、测试和预构建运行时。
- offline/runtime.tar.gz：仅预构建 Conda 环境。工作站有网络和 Conda 时可以完全不复制它，直接执行 bash setup.sh，脚本会根据 environment.yml 重建环境。

复制大型文件到工作站后，先运行：
  sha256sum -c large_files/portable/lesion_complete.tar.gz.sha256
  sha256sum -c large_files/offline/runtime.tar.gz.sha256

setup.sh 优先查找 large_files/offline/runtime.tar.gz，其次兼容旧路径 offline/runtime.tar.gz。GitHub 克隆后没有大型快照也可以通过 environment.yml 安装依赖。
