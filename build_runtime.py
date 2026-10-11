#!/usr/bin/env python
"""构建可移植 Conda 运行环境；修正 conda-pack 重写硬链接成员时的 tar 元数据。

修正只作用于本次打包进程，不修改 env 内安装包。先完整验证，再原子提交新归档。
"""
from __future__ import annotations

import argparse
from io import BytesIO
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time
import conda_pack
from conda_pack.formats import TarArchive
from functions.storage import digest, save_json


def add_regular_bytes(self, source: str, sourcebytes: bytes, target: str) -> None:
    """重写内容必须用普通文件头；硬链接头不能携带数据体。"""
    info = self.archive.gettarinfo(source, target)
    # 同 inode 已收录时 gettarinfo 可能返回 LNKTYPE；清除链接语义后才写数据。
    info.type = tarfile.REGTYPE
    info.linkname = ""
    info.pax_headers.pop("linkpath", None)
    info.size = len(sourcebytes)
    self.archive.addfile(info, BytesIO(sourcebytes))


def hardlink_regression(root: Path) -> dict:
    """最小回归：两条路径指向同一 inode，第二个内容经 add_bytes 重写。"""
    directory = root / "archive/upgrade"
    directory.mkdir(parents=True, exist_ok=True)
    buffer = BytesIO()
    # 自建硬链接夹具，使首次环境及移植后已解包环境都能重复执行回归。
    with tempfile.TemporaryDirectory(prefix="hardlink_probe_", dir=directory) as temporary:
        source = Path(temporary) / "first"
        partner = Path(temporary) / "second"
        content = b"synthetic activation script\n"
        source.write_bytes(content)
        os.link(source, partner)
        source_inode = source.stat().st_ino
        with TarArchive(buffer, "") as archive:
            archive.add(str(source), "first.fish")
            archive.add_bytes(str(partner), content, "second.fish")
            archive.add_bytes(str(partner), b"sentinel\n", "last.txt")
    buffer.seek(0)
    with tarfile.open(fileobj=buffer, mode="r:") as archive:
        members = archive.getmembers()
        if [item.name for item in members] != ["first.fish", "second.fish", "last.txt"]:
            raise ValueError("硬链接归档导致后续成员丢失")
        if not all(item.isfile() for item in members):
            raise ValueError("带数据的成员仍为链接头")
        if archive.extractfile("second.fish").read() != content or archive.extractfile("last.txt").read() != b"sentinel\n":
            raise ValueError("最小硬链接回归内容不一致")
    return {"status":"passed", "members":3, "source_inode":source_inode,
            "corrected_member_type":"REGTYPE", "environment_files_modified":False}


def verify_tar(path: Path, log_dir: Path) -> dict:
    """GNU tar 和 Python 独立遍历，成员数必须一致；仅通过 gzip 校验不够。"""
    log_dir.mkdir(parents=True, exist_ok=True)
    names_path = log_dir / "members.txt"
    errors_path = log_dir / "tar_errors.txt"
    with names_path.open("w") as names, errors_path.open("w") as errors:
        result = subprocess.run(["tar", "-tzf", str(path)], stdout=names, stderr=errors, check=False)
    if result.returncode != 0 or errors_path.stat().st_size:
        raise ValueError(f"GNU tar 检查失败：{errors_path}")
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        invalid = [item.name for item in members if (item.islnk() or item.issym()) and item.size]
        if invalid:
            raise ValueError(f"链接成员带非零数据长度：{invalid[:5]}")
    listed = len(names_path.read_text().splitlines())
    if listed != len(members):
        raise ValueError("GNU tar 与 Python 读取的归档成员数不一致")
    return {"gnu_tar_returncode":0, "gnu_tar_stderr_bytes":0, "members":len(members),
            "last_member":members[-1].name, "nonzero_link_members":0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--replace-invalid", action="store_true", help="确认旧包无效后保留到 archive/upgrade/offline_invalid")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    offline = root / "offline"
    output = offline / "runtime.tar.gz"
    temporary = offline / "runtime.tar.gz.part"
    if temporary.exists():
        raise FileExistsError("上一次 .part 已存在，请核查后使用独立新位置")
    if args.threads < 1 or args.threads > 32:
        raise ValueError("threads 必须位于 1–32")
    # 一、保留已确认无效的旧包，不触碰环境或历史原版运行包。
    if output.exists():
        if not args.replace_invalid:
            raise FileExistsError("运行包已存在；不会覆盖有效离线包")
        probe = subprocess.run(["tar", "-tzf", str(output)], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=False)
        if probe.returncode == 0:
            raise ValueError("现有运行包通过 GNU tar 检查，拒绝按无效包替换")
        backup = root / "archive/upgrade/offline_invalid"
        backup.mkdir(parents=True, exist_ok=False)
        backup.resolve().relative_to(root.resolve())
        old_hash = digest(output)
        output.rename(backup / output.name)
        checksum = output.with_name(output.name + ".sha256")
        if checksum.exists():
            checksum.rename(backup / checksum.name)
        save_json(backup / "diagnosis.json", {"sha256":old_hash,"reason":"hardlink header with nonzero size and body in conda_pack TarArchive._add_bytes",
                  "gnu_tar_returncode":probe.returncode,"gnu_tar_stderr":probe.stderr.decode(errors="replace")})
    # 二、当前进程局部修正；不执行新子进程，从而保留 monkeypatch 的作用。
    original = TarArchive._add_bytes
    started = time.monotonic()
    try:
        TarArchive._add_bytes = add_regular_bytes
        regression = hardlink_regression(root)
        print(json.dumps({"hardlink_regression":regression}), flush=True)
        conda_pack.pack(prefix=str(root / "env"), output=str(temporary), format="tar.gz",
                        n_threads=args.threads, compress_level=4, verbose=True, force=False)
    finally:
        TarArchive._add_bytes = original
    # 三、双读取器验证通过后提交；新 SHA 与构建证据同时保留。
    verification = verify_tar(temporary, root / "archive/upgrade/runtime_validation")
    temporary.rename(output)
    checksum = digest(output)
    output.with_name(output.name + ".sha256").write_text(f"{checksum}  {output.name}\n")
    record = {"status":"passed","conda_pack_version":conda_pack.__version__,
              "fix":"TarArchive._add_bytes writes REGTYPE with empty linkname", "environment_packages_modified":False,
              "threads":args.threads,"compression_level":4,"elapsed_seconds":time.monotonic()-started,
              "bytes":output.stat().st_size,"sha256":checksum,"regression":regression,"tar_verification":verification}
    save_json(offline / "build.json", record)
    print(json.dumps(record,ensure_ascii=False),flush=True)


if __name__ == "__main__":
    main()
