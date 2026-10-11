#!/usr/bin/env bash
# 项目环境安装入口：优先使用可选离线 conda-pack；没有离线包时可由 Conda 根据 environment.yml 在线重建。
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" # 项目绝对路径
runtime_candidates=("$project_dir/large_files/offline/runtime.tar.gz" "$project_dir/offline/runtime.tar.gz")
runtime_path=""
for candidate in "${runtime_candidates[@]}"; do
  if [[ -f "$candidate" ]]; then runtime_path="$candidate"; break; fi
done
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo '仅支持 Linux x86_64；本交付在 Ubuntu 20.04 验证。' >&2
  exit 2
fi
glibc_version="$(getconf GNU_LIBC_VERSION | cut -d' ' -f2)" # 目标系统 glibc 版本
if [[ "$(printf '%s\n' 2.31 "$glibc_version" | sort -V | head -n1)" != 2.31 ]]; then
  echo '需要 glibc >= 2.31。' >&2
  exit 2
fi
if [[ -f env/.runtime_ready ]]; then
  if [[ "$(cat env/.runtime_ready)" != "$project_dir" ]]; then
    echo '已安装环境被直接移动；请在新位置重新解包或重建。' >&2
    exit 2
  fi
elif [[ -e env ]]; then
  echo 'env 已存在但没有本项目安装标识；保留现有目录，请先检查。' >&2
  exit 2
elif [[ -n "$runtime_path" ]]; then
  checksum_path="${runtime_path}.sha256"
  [[ -f "$checksum_path" ]] || { echo "缺少离线包校验文件: $checksum_path" >&2; exit 2; }
  (cd "$(dirname "$runtime_path")" && sha256sum -c "$(basename "$checksum_path")")
  mkdir "$project_dir/env"
  tar --no-same-owner -xzf "$runtime_path" -C "$project_dir/env"
  "$project_dir/env/bin/python" "$project_dir/env/bin/conda-unpack"
  printf '%s\n' "$project_dir" > "$project_dir/env/.runtime_ready"
else
  # 没有预构建运行时，使用工作站 Conda 按 environment.yml 建立同名项目环境。
  conda_command="${CONDA_EXE:-}"
  if [[ -z "$conda_command" ]]; then conda_command="$(command -v conda || true)"; fi
  [[ -n "$conda_command" ]] || { echo '未找到 Conda；请提供 large_files/offline/runtime.tar.gz 或安装 Conda。' >&2; exit 2; }
  "$conda_command" env create --prefix "$project_dir/env" --file "$project_dir/environment.yml" --yes
  printf '%s\n' "$project_dir" > "$project_dir/env/.runtime_ready"
fi
"$project_dir/env/bin/python" -B -c 'import numpy, scipy, nibabel, h5py, numba, pandas, matplotlib, yaml, psutil, pyreadstat; print("Python 环境可用")'
"$project_dir/env/bin/Rscript" --vanilla -e 'for (package in c("ordinal","mice","mitml","ordinalNet","sandwich","lme4","jsonlite")) stopifnot(requireNamespace(package, quietly=TRUE)); cat("R 环境可用\n")'
