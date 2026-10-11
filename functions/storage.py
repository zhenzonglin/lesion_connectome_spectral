"""一、文件校验；二、原子提交；三、项目内路径约束。"""

from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path


def digest(path: Path) -> str:
    """分块计算文件 SHA256，不将大型索引一次读入内存。"""
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def fingerprint(value: object) -> str:
    """对参数/文件清单生成与字典顺序无关的指纹。"""
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, ensure_ascii=False).encode()
    ).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, value: dict) -> None:
    """先写同目录临时文件并同步，再替换；不用于人工审核文件。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def inside(root: Path, relative: str) -> Path:
    """用户路径必须位于项目内，禁止 ../ 和外部软链接。"""
    path = (root / relative).resolve()
    path.relative_to(root.resolve())
    return path


def verify_files(root: Path, records: dict) -> None:
    for name, expected in records.items():
        path = inside(root, name)
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"文件缺失或校验改变：{path}")
