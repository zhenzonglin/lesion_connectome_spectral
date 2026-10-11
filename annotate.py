#!/usr/bin/env python
"""按需导出三种空间对象或运行已验证组合的空间注释；不改变临床模型。"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from functions.configuration import load_config
from functions.annotation import export_healthy_modes, export_node_map, loss_contributions, run_annotation
from functions.reference import load_reference, verify_reference
from functions.storage import save_json, read_json, verify_files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--output", required=True, help="项目内尚不存在的输出目录")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--healthy", action="store_true", help="导出所有129个健康模态")
    selection.add_argument("--patient", help="某例紧凑result.npz路径，必须是有效计算结果")
    selection.add_argument("--input", help="包含maps_manifest.json的已验证组合脑图目录")
    parser.add_argument("--mode", type=int, help="患者按需导出的模态编号，1至129")
    parser.add_argument("--nifti", action="store_true", help="额外保存NIfTI；默认只有PNG和脑区表")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    config = load_config(root / args.config)
    output = (root / args.output).resolve()
    output.relative_to(root)
    if output.exists():
        raise FileExistsError("导出拒绝覆盖既有目录，请使用新名称")
    if args.input:
        input_dir = (root / args.input).resolve()
        input_dir.relative_to(root)
        print(run_annotation(root, config, input_dir, output))
        return
    manifest = verify_reference(root)
    reference = load_reference(root, with_index=False)
    figure = config.get("annotation", {}).get("figure", {})
    if args.healthy:
        print(export_healthy_modes(reference, output, args.nifti, figure))
        return
    if args.mode not in reference["mode_ids"]:
        raise ValueError("患者空间图需要 --mode 1至129")
    source = (root / args.patient).resolve()
    source.relative_to(root)
    record = read_json(source.parent / "record.json")
    if record.get("automatic_status") != "passed" or record.get("reference_id") != manifest["reference_id"]:
        raise ValueError("患者结果未通过计算验收或引用了不同健康参考")
    verify_files(source.parent, record["outputs"])
    with np.load(source) as patient:
        local, network = loss_contributions(reference, patient["lost_counts"], patient["roi_coverage"])
    position = int(np.flatnonzero(reference["mode_ids"] == args.mode)[0])
    for channel, matrix in (("g", local), ("w", network)):
        export_node_map(reference, matrix[:, position], output / channel,
                        f"Mode {args.mode}: {channel.upper()} loss contribution", args.nifti, figure)
    save_json(output / "description.json", {"mode": args.mode,
              "source": str(source.relative_to(root)), "g_sum": float(np.nansum(local[:, position])),
              "w_sum": float(np.nansum(network[:, position])),
              "interpretation": "节点损伤贡献；不是患者重算的特征向量"})
    print(output)


if __name__ == "__main__":
    main()
