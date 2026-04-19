"""Draw DAGs for base/quant open-economy HANK models.

Usage examples:
  python draw_hank_dag.py --model base
  python draw_hank_dag.py --model quant
"""

from __future__ import annotations

import argparse
from pathlib import Path

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

    print(f"DOT: {out_dot}")
    if used_graphviz:
        print(f"PNG (drawdag/graphviz): {out_base.with_suffix('.png')}")

if __name__ == "__main__":
    main()
