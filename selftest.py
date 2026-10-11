#!/usr/bin/env python
"""移植后自检：验证参考校验值、数值定义、完整线段算法与模拟病灶。"""
import os
import sys
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"]:
        os.environ[key] = "1"
    os.environ["NUMBA_CACHE_DIR"] = str(root / ".cache/numba")
    from functions.reference import verify_reference, numerical_checks, load_reference

    print(numerical_checks(load_reference(root, with_index=False)), flush=True)
    verify_reference(root)
    import pytest
    code = pytest.main(
        [
            str(root / "tests"),
            "-q",
            "-p",
            "no:cacheprovider",
        ]
    )
    if code: return int(code)
    import subprocess
    if '--statistics' in sys.argv:
        output=root/'archive/offline_selftest'
        if output.exists(): raise FileExistsError('统计自检目录已存在，请保留后使用新的解包项目运行')
        return subprocess.run([str(root/'env/bin/Rscript'),'--vanilla',str(root/'statistics/selftest.R'),str(output)]).returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
