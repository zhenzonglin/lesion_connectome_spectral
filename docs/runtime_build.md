# 离线环境构建与归档格式修正

`build_runtime.py` 使用当前项目环境内的 conda-pack API，在同一个Python进程中构建 `offline/runtime.tar.gz`。它不改动安装包源码，不删除环境，不依赖恢复原始包缓存。目标机器仍由 `setup.sh` 解包并执行 conda-unpack。

本轮发现一种可复现的归档问题：`env/bin/activate.fish` 与 conda-pack 包内的同名脚本是硬链接。`TarArchive._add_bytes` 从Python `gettarinfo()`取得硬链接类型后改写 `size` 并写入数据，却没有把类型改回普通文件。tar规则中硬链接没有文件数据体。GNU tar因此报告跳到下一个头；Python在脚本数据处停止读取。gzip校验和文件SHA256均可正常通过，不能代替tar结构检查。

修正仅在本次打包进程中替换 `TarArchive._add_bytes`：保留权限和时间等元数据，将携带实际数据的条目标记为 `REGTYPE`，清空 `linkname/linkpath` 后写入字节。未改写内容的正常硬链接仍按原方式保存；修正不是压缩线程设置。

构建顺序：

1. 若指定 `--replace-invalid`，先用GNU tar确认旧包无效，再将旧包及SHA移入 `archive/upgrade/offline_invalid`，保存原因和校验值；不覆盖健康历史包。
2. 运行三成员硬链接最小回归，检查重写内容及最后成员仍可读取。
3. conda-pack打包到 `.part`，默认4线程、压缩级别4。
4. GNU tar无错误遍历；Python独立完整遍历；两者成员数相同；所有链接成员数据长度为0。
5. 原子提交新包及SHA256，写入 `offline/build.json`。不同路径离线解包与Python/R执行由移植验收另行验证。

```bash
# 初次构建（offline/runtime.tar.gz不存在时）
./env/bin/python build_runtime.py

# 只允许替换已被GNU tar判定无效的包；旧包原样保留
./env/bin/python build_runtime.py --replace-invalid
```

独立回归代码为 `tests/test_runtime_pack.py`。构建过程使用源码可查的本机 conda-pack，实际版本号与修正说明在 `offline/build.json`；完整归档成员表与GNU tar诊断保存在 `archive/upgrade/runtime_validation`。此归档格式修正不改变健康底座、患者特征或统计结果。
