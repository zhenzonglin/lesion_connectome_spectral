"""离线打包回归：硬链接经内容改写后必须成为普通成员，后续成员仍能读取。"""
from io import BytesIO
import os
import tarfile
from conda_pack.formats import TarArchive
from build_runtime import add_regular_bytes, hardlink_regression


def test_rewritten_hardlink_keeps_archive_readable(tmp_path, monkeypatch):
    source = tmp_path / "first"
    partner = tmp_path / "linked"
    source.write_bytes(b"original bytes\n")
    os.link(source, partner)
    monkeypatch.setattr(TarArchive, "_add_bytes", add_regular_bytes)
    buffer = BytesIO()
    with TarArchive(buffer, "") as archive:
        archive.add(str(source), "first")
        archive.add_bytes(str(partner), b"modified independent bytes\n", "rewritten")
        archive.add_bytes(str(partner), b"tail\n", "last")
        archive.add(str(partner), "plain_hardlink")
    buffer.seek(0)
    with tarfile.open(fileobj=buffer, mode="r:") as archive:
        assert archive.getnames() == ["first","rewritten","last","plain_hardlink"]
        assert archive.getmember("rewritten").isreg()
        assert archive.extractfile("rewritten").read() == b"modified independent bytes\n"
        assert archive.extractfile("last").read() == b"tail\n"
        assert archive.getmember("plain_hardlink").islnk()
        assert archive.extractfile("plain_hardlink").read() == source.read_bytes()


def test_builder_regression_does_not_need_original_environment_links(tmp_path, monkeypatch):
    monkeypatch.setattr(TarArchive, "_add_bytes", add_regular_bytes)
    record = hardlink_regression(tmp_path)
    assert record["status"] == "passed"
    assert record["members"] == 3
    assert not (tmp_path / "env").exists()
