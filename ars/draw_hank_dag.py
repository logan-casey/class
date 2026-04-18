"""Draw DAGs for base/quant open-economy HANK models.

Usage examples:
  python draw_hank_dag.py --model base
  python draw_hank_dag.py --model quant

If `sequence_jacobian.drawdag` + Graphviz are available, this script uses them.
Otherwise, it falls back to a matplotlib-based DAG figure and always writes a DOT file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import arss_hank_base as base_mod
import arss_hank_quant as quant_mod


BASE_SPEC = {
    "inputs": ["i_star"],
    "unknowns": ["Y", "W", "Q", "E", "i", "r"],
    "targets": ["goods_mkt", "wnkpc", "exrate_res", "mp_res", "fisher_res", "uip_res"],
}

QUANT_SPEC = {
    "inputs": ["i_star"],
    "unknowns": ["Y", "W", "PH", "PH_star", "Q", "E", "i", "r", "xH_hat", "xH_star_hat", "shareH", "shareH_star"],
    "targets": [
        "goods_mkt",
        "wnkpc",
        "piH_res",
        "piH_star_res",
        "mp_taylor_res",
        "fisher_res",
        "uip_res",
        "xH_target_res",
        "xH_star_target_res",
        "shareH_res",
        "shareH_star_res",
        "exrate_res",
    ],
}


def get_model_and_spec(name: str):
    if name == "base":
        model, _ = base_mod.build_base_model()
        return model, BASE_SPEC
    if name == "quant":
        model, _ = quant_mod.build_quant_model()
        return model, QUANT_SPEC
    raise ValueError(f"Unknown model '{name}'. Use 'base' or 'quant'.")


def write_dot(model, inputs, unknowns, targets, out_dot: Path) -> None:
    lines = [
        "digraph HANK_DAG {",
        "  rankdir=LR;",
        '  exog [shape=box,label="exogenous"];',
        '  unk [shape=box,label="unknowns"];',
        '  tar [shape=diamond,label="targets"];',
    ]
    for i, b in enumerate(model.blocks):
        lines.append(f'  b{i} [label="{i}: {b.name}"];')

    for i, b in enumerate(model.blocks):
        b_inputs = set(b.inputs)
        from_exog = b_inputs & set(inputs)
        from_unk = b_inputs & set(unknowns)
        if from_exog:
            lines.append(f'  exog -> b{i} [label="{", ".join(sorted(from_exog))}"];')
        if from_unk:
            lines.append(f'  unk -> b{i} [label="{", ".join(sorted(from_unk))}"];')

        for j in model.revadj[i]:
            lab = b_inputs & set(model.blocks[j].outputs)
            lines.append(f'  b{j} -> b{i} [label="{", ".join(sorted(lab))}"];')

        for t in targets:
            if t in b.outputs:
                lines.append(f'  b{i} -> tar [label="{t}"];')

    lines.append("}")
    out_dot.write_text("\n".join(lines), encoding="utf-8")


def draw_custom(model, inputs, unknowns, targets, out_png: Path) -> None:
    n = len(model.blocks)
    x_exog, x_unk, x_blk, x_tar = 0.05, 0.22, 0.52, 0.9
    y_top, y_bot = 0.95, 0.05
    ys = [y_top - i * (y_top - y_bot) / max(1, n - 1) for i in range(n)]

    fig, ax = plt.subplots(figsize=(15, 11))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def draw_box(x, y, text, w=0.22, h=0.04, fc="#F5F7FA", ec="#4B5563", fs=8):
        rect = FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.01",
            linewidth=0.8,
            edgecolor=ec,
            facecolor=fc,
        )
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=fs)

    def draw_arrow(x1, y1, x2, y2, label="", color="#6B7280", alpha=0.7):
        arr = FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="->",
            mutation_scale=8,
            linewidth=0.8,
            color=color,
            alpha=alpha,
            connectionstyle="arc3,rad=0",
        )
        ax.add_patch(arr)
        if label:
            xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(xm, ym + 0.008, label, fontsize=6, color=color, ha="center", va="bottom")

    draw_box(x_exog, 0.90, "exogenous")
    draw_box(x_unk, 0.90, "unknowns")
    draw_box(x_tar, 0.90, "targets", fc="#FEF3C7", ec="#92400E")

    for i, b in enumerate(model.blocks):
        btype = "HA" if "het" in b.__class__.__name__.lower() else "block"
        draw_box(x_blk, ys[i], f"{i}: {b.name} [{btype}]", w=0.34, h=0.035, fs=7)

    for i, b in enumerate(model.blocks):
        b_inputs = set(b.inputs)
        ex_l = b_inputs & set(inputs)
        un_l = b_inputs & set(unknowns)
        if ex_l:
            draw_arrow(x_exog + 0.11, 0.90, x_blk - 0.18, ys[i], ",".join(sorted(ex_l)), color="#2563EB", alpha=0.6)
        if un_l:
            draw_arrow(x_unk + 0.11, 0.90, x_blk - 0.18, ys[i], ",".join(sorted(un_l)), color="#7C3AED", alpha=0.6)

        for j in model.revadj[i]:
            lab = b_inputs & set(model.blocks[j].outputs)
            draw_arrow(x_blk + 0.17, ys[j], x_blk - 0.17, ys[i], ",".join(sorted(lab)), color="#6B7280", alpha=0.4)

        for t in targets:
            if t in b.outputs:
                draw_arrow(x_blk + 0.18, ys[i], x_tar - 0.11, 0.90, t, color="#B45309", alpha=0.7)

    ax.set_title(f"{model.name} DAG (custom fallback)", fontsize=12)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def try_graphviz_drawdag(model, inputs, unknowns, targets, out_base_no_suffix: Path) -> bool:
    try:
        from sequence_jacobian import drawdag
    except Exception:
        return False

    try:
        drawdag(model, inputs, unknowns, targets, leftright=True, save=True, savepath=str(out_base_no_suffix))
        return out_base_no_suffix.with_suffix(".png").exists()
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["base", "quant"], default="base")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    model, spec = get_model_and_spec(args.model)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_base = outdir / f"ha_oe_{args.model}_dag"
    out_dot = out_base.with_suffix(".dot")
    out_png_custom = outdir / f"ha_oe_{args.model}_dag_custom.png"

    write_dot(model, spec["inputs"], spec["unknowns"], spec["targets"], out_dot)
    used_graphviz = try_graphviz_drawdag(model, spec["inputs"], spec["unknowns"], spec["targets"], out_base)
    if not used_graphviz:
        draw_custom(model, spec["inputs"], spec["unknowns"], spec["targets"], out_png_custom)

    print(f"DOT: {out_dot}")
    if used_graphviz:
        print(f"PNG (drawdag/graphviz): {out_base.with_suffix('.png')}")
    else:
        print(f"PNG (custom fallback): {out_png_custom}")


if __name__ == "__main__":
    main()
