#!/usr/bin/env python
"""按需导出健康模式、患者 SC 或模式损失贡献；拒绝覆盖同名输出。"""
import argparse
from pathlib import Path
import nibabel as nib
import numpy as np
from functions.reference import verify_reference, load_reference
from functions.spectral import fixed_support_loss,tissue_loss
from functions.plotting import mode_png
from functions.storage import save_json, read_json, fingerprint


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", type=int)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--channel", choices=["g","w"], default="w")
    parser.add_argument("--squared", action="store_true", help="导出健康平方载荷")
    parser.add_argument(
        "--mapped-lesion", action="store_true", help="仅导出统一空间病灶供人工检查"
    )
    parser.add_argument("--healthy", action="store_true")
    parser.add_argument("--run")
    parser.add_argument("--subject")
    parser.add_argument("--name", help="新的导出目录短名，默认按对象命名")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    from functions.configuration import load_config
    from functions.control import run_directory
    config = load_config(root / args.config)
    base = config["paths"]["results"]
    def patient_directory():
        directory = run_directory(root,args.run,base)
        if (directory/'stages.json').exists():
            state = read_json(directory/'stages.json')
            if state['stages'].get('lesion',{}).get('status') != 'completed':
                raise ValueError('病灶阶段未完成')
            directory=Path(state['stages']['lesion']['path'])
        target=(directory/'subjects'/args.subject).resolve()
        target.relative_to((directory/'subjects').resolve())
        return target
    manifest = verify_reference(root)
    reference = load_reference(root, with_index=False)
    if args.mapped_lesion:
        if not args.run or not args.subject or args.healthy:
            raise ValueError("病灶导出需要 --run 和 --subject")
        from functions.control import run_directory
        from functions.storage import digest, inside
        from functions.space import map_patient

        source_directory = patient_directory()
        record = read_json(source_directory / "record.json")
        source = root / record["input_path"]
        if digest(source) != record["input_sha256"]:
            raise ValueError("原始病灶已改变，不能重建旧结果")
        if record["reference_id"] != manifest["reference_id"]:
            raise ValueError("健康参考版本不同")
        mask, summary = map_patient(source, reference["template"])
        output = run_directory(
            root, args.name or (args.subject + "_mapped_lesion"), base
        )
        output.mkdir(parents=True, exist_ok=False)
        nib.save(
            nib.Nifti1Image(mask.astype(np.uint8), reference["template"].affine),
            output / "lesion.nii.gz",
        )
        save_json(output / "description.json", summary)
        print(output)
        return
    if args.mode is None:
        raise ValueError("模式导出需要 --mode；病灶导出使用 --mapped-lesion")
    selected = np.flatnonzero(reference["mode_ids"] == args.mode)
    if len(selected) != 1:
        raise ValueError("请选择 1—129 范围内的非零模态")
    position = int(selected[0])
    if args.healthy:
        values = reference["vectors"][:, position].copy()
        if args.squared: values = values**2
        label = f"healthy_mode_{args.mode}"
        kind = "健康模式有符号空间载荷"
    else:
        if not args.run or not args.subject:
            raise ValueError("患者导出需要 --run 和 --subject")
        from functions.control import run_directory

        target = patient_directory()
        record = read_json(target / "record.json")
        if record["automatic_status"] != "passed":
            raise ValueError("该患者已排除，没有可导出的模式特征")
        if record["reference_id"] != manifest["reference_id"]:
            raise ValueError("患者与当前健康参考不同")
        from functions.storage import verify_files

        verify_files(target, record["outputs"])
        with np.load(target / "result.npz") as data:
            lost = data["lost_counts"].copy()
            rho = data["roi_coverage"].copy()
        ratios, _, spatial = fixed_support_loss(lost, reference, True)
        if args.channel == 'g':
            ratios, spatial = tissue_loss(rho,reference,True)
        values = spatial[:, position]
        label = f"{args.subject}_{args.channel}_mode_{args.mode}"
        kind = "该模态的脑区损失贡献；脑区值求和等于支撑损失"
    from functions.control import run_directory

    output = run_directory(root, args.name or label, base)
    output.mkdir(parents=True, exist_ok=False)
    values[~reference["active"]] = np.nan
    mapped = np.r_[np.nan, values][reference["atlas"]]
    nib.save(
        nib.Nifti1Image(mapped.astype(np.float32), reference["template"].affine),
        output / "map.nii.gz",
    )
    table = reference["nodes"][["index", "label", "group"]].copy()
    table["value"] = values
    table.to_csv(output / "regions.tsv", sep="\t", index=False)
    mode_png(
        reference["template"], reference["atlas"], values, label, output / "map.png",
        dpi=config['annotation']['figure']['dpi'],
        cmap_signed=config['annotation']['figure']['colormap_signed'],
        cmap_loss=config['annotation']['figure']['colormap_loss']
    )
    if not args.healthy:
        graph = reference["weights"].copy()
        i, j = reference["edge_i"], reference["edge_j"]
        remaining = (reference["counts"][i, j] - lost) / (
            reference["volumes"][i] + reference["volumes"][j]
        )
        graph[i, j] = remaining
        graph[j, i] = remaining
        np.savez_compressed(output / "sc.npz", weights=graph)
    save_json(
        output / "description.json",
        {
            "object": kind,
            "mode": args.mode,
            "excluded_nodes": "NaN",
            "voxel_values": "constant parcel value; not voxel density; do not sum voxels as node contributions",
        },
    )
    print(output)


if __name__ == "__main__":
    main()
