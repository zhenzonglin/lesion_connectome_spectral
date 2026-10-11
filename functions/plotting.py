"""必要 PNG 质控和按需脑图；使用固定方向标记，避免左右混淆。"""

from __future__ import annotations
from pathlib import Path
import nibabel as nib
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cohort_spectrum(source: Path, path: Path) -> None:
    """小样本逐例曲线；全量以中位数和四分位带显示，不绘制数万条曲线。"""
    import pandas as pd

    table = pd.read_csv(source, sep="\t")
    columns = [name for name in table.columns if name.startswith('w_') or name.startswith('mode_')]
    values = table[columns].to_numpy()
    modes = np.array([int(name.split("_")[1]) for name in columns])
    figure, axis = plt.subplots(figsize=(11, 4), dpi=120)
    if len(table) <= 12:
        for subject, row in zip(table.subject_id, values):
            axis.plot(modes, row, label=subject, lw=1.2)
    else:
        low, middle, high = np.quantile(values, [0.25, 0.5, 0.75], axis=0)
        axis.fill_between(modes, low, high, alpha=0.2, label="Interquartile range")
        axis.plot(modes, middle, label=f"Median (n={len(table)})")
    axis.set(
        xlabel="Healthy mode index",
        ylabel="Structural support loss (fraction)",
        ylim=(0, None),
    )
    axis.grid(alpha=0.2)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def reference_png(reference: dict, path: Path) -> None:
    """参考概览：健康权重、强度、完整非零谱和节点范围。"""
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=120)
    weights = np.ma.masked_where(reference["weights"] == 0, reference["weights"])
    image = axes[0, 0].imshow(np.ma.log10(weights), cmap="viridis")
    figure.colorbar(image, ax=axes[0, 0], label="log10(weight)")
    axes[0, 0].set(
        title="Frozen weighted SC",
        xlabel="ROI index (zero based)",
        ylabel="ROI index (zero based)",
    )
    axes[0, 1].plot(np.arange(1, 136), reference["strength"], lw=1)
    axes[0, 1].set(
        title="Healthy node strength", xlabel="ROI label", ylabel="Weighted strength"
    )
    axes[1, 0].plot(reference["mode_ids"], reference["values"])
    axes[1, 0].set(
        title="All nonzero healthy modes", xlabel="Mode index", ylabel="Eigenvalue"
    )
    excluded = reference["nodes"].loc[~reference["active"], "label"].to_list()
    axes[1, 1].axis("off")
    axes[1, 1].text(
        0,
        0.95,
        "LQT135: 135 anatomical / 130 connected nodes\n129 nonzero modes / 1375 edges\n\nExcluded from harmonics:\n"
        + "\n".join(excluded)
        + "\n\nWhole brainstem = mixed tissue\nAnatomical review: PENDING",
        va="top",
        fontsize=11,
    )
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def mode_png(
    template,
    atlas,
    values,
    title: str,
    path: Path,
    note: str = "parcel loadings; excluded nodes masked",
    dpi: int = 180,
    cmap_signed: str = "coolwarm",
    cmap_loss: str = "magma",
) -> None:
    """按需将脑区数值映回空间；排除节点为 NaN，不人为插值。"""
    image = nib.as_closest_canonical(
        nib.Nifti1Image(
            np.r_[np.nan, values][atlas].astype(np.float32), template.affine
        )
    )
    background = np.asarray(nib.as_closest_canonical(template).dataobj)
    data = np.asarray(image.dataobj)
    bound = float(np.nanmax(np.abs(values))) or 1.0
    nonnegative = bool(np.nanmin(values) >= 0)
    figure, axes = plt.subplots(2, 4, figsize=(14, 7), dpi=dpi)
    for axis, z_mm in zip(axes.ravel(), [-45, -25, -5, 15, 35, 55, 65, 75]):
        index = int(round((z_mm - image.affine[2, 3]) / image.affine[2, 2]))
        axis.imshow(background[:, :, index].T, cmap="gray", origin="lower")
        artist = axis.imshow(
            np.ma.masked_invalid(data[:, :, index].T),
            cmap=cmap_loss if nonnegative else cmap_signed,
            origin="lower",
            vmin=0 if nonnegative else -bound,
            vmax=bound,
            alpha=0.9,
        )
        axis.set_title(f"z={z_mm} mm | L / R")
        axis.axis("off")
    figure.suptitle(title + " | " + note)
    figure.subplots_adjust(right=0.9, top=0.88)
    color_axis = figure.add_axes([0.92, 0.2, 0.015, 0.6])
    figure.colorbar(artist, cax=color_axis, label="Parcel value")
    figure.savefig(path)
    plt.close(figure)
