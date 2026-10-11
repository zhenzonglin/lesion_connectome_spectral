#!/usr/bin/env python
"""最小离线移植包：白名单收集、路径/链接检查、哈希验证、同名不覆盖。"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import tarfile
from pathlib import Path
from functions.storage import digest, save_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", default="lesion_spectral")
    parser.add_argument("--include-healthy-images", action="store_true", help="连同可由数值基底再生的258张健康图打包")
    args = parser.parse_args()
    if not args.name.replace("_", "").isalnum():
        raise ValueError("包名只允许字母、数字和下划线")
    root = Path(__file__).resolve().parent
    output = root / "portable"
    output.mkdir(exist_ok=True)
    archive = output / (args.name + ".tar.gz")
    if archive.exists() or archive.with_suffix(".part").exists():
        raise FileExistsError("移植包已存在，指定新 --name")
    paths = (
        list(root.glob("*.py"))
        + list(root.glob("*.sh"))
        + list(root.glob("*.md"))
        + list(root.glob("*.yml"))
        + list(root.glob("*.yaml"))
    )
    for folder in ["functions", "reference", "tests", "docs", "statistics"]:
        paths += [
            p
            for p in (root / folder).rglob("*")
            if p.is_file()
            and "__pycache__" not in p.parts
            and ".pytest_cache" not in p.parts
        ]
    runtime_root = root / "large_files/offline"
    if not (runtime_root / "runtime.tar.gz").is_file():
        runtime_root = root / "offline"
    paths += [
        runtime_root / "runtime.tar.gz",
        runtime_root / "runtime.tar.gz.sha256",
        root / "offline/versions.json",
        root / "offline/build.json",
        root / "clinical/readme.md",
        root / "qc/code_checklist.tsv",
        root / "qc/readme.md",
        root / "qc/reference.png",
        root / "qc/anatomy.png",
    ]
    records = {}
    for path in sorted(set(paths)):
        if not args.include_healthy_images and path.is_relative_to(root/'reference/annotation/healthy'):
            continue
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.is_symlink():
            raise ValueError(f"工程资源不能依赖软链接：{path}")
        relative = str(path.relative_to(root))
        records[relative] = {"sha256": digest(path), "bytes": path.stat().st_size}
    manifest = {
        "files": records,
        "excluded": [
            "archive",
            "clinical/*.sas7bdat",
            "clinical/ID mappings and member lists",
            "lesion",
            "results",
            "env",
            "checks",
            "qc/manual_review.tsv",
            ".cache",
        ],
        "platform": "Linux x86_64; Ubuntu >=20.04; glibc >=2.31",
        "patient_data_included": False,
        "healthy_images_included": args.include_healthy_images,
        "healthy_images_regeneration": "env/bin/python annotate.py --healthy --output reference/annotation/healthy",
    }
    base = "lesion_connectome_spectral"
    with tarfile.open(archive.with_suffix(".part"), "w:gz", compresslevel=6) as stream:
        for relative in records:
            stream.add(root / relative, arcname=f"{base}/{relative}", recursive=False)
        for folder in ["lesion", "results", "qc", "clinical"]:
            entry = tarfile.TarInfo(f"{base}/{folder}")
            entry.type = tarfile.DIRTYPE
            entry.mode = 0o755
            stream.addfile(entry)
        content = json.dumps(manifest, ensure_ascii=False, indent=2).encode()
        entry = tarfile.TarInfo(f"{base}/package_manifest.json")
        entry.size = len(content)
        entry.mode = 0o644
        stream.addfile(entry, io.BytesIO(content))
    archive.with_suffix(".part").replace(archive)
    with tarfile.open(archive, "r:gz") as stream:
        for member in stream:
            if member.isfile() and member.name != f"{base}/package_manifest.json":
                relative = member.name[len(base) + 1 :]
                value = hashlib.sha256(stream.extractfile(member).read()).hexdigest()
                if value != records[relative]["sha256"]:
                    raise ValueError(f"归档校验失败：{relative}")
    checksum = digest(archive)
    archive.with_name(archive.name + ".sha256").write_text(
        f"{checksum}  {archive.name}\n"
    )
    save_json(
        output / (args.name + "_manifest.json"),
        {
            **manifest,
            "archive_sha256": checksum,
            "archive_bytes": archive.stat().st_size,
        },
    )
    print(
        json.dumps(
            {
                "archive": str(archive),
                "bytes": archive.stat().st_size,
                "sha256": checksum,
                "files": len(records),
            }
        )
    )


if __name__ == "__main__":
    main()
