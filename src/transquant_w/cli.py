"""Command-line interface: transquant-w --target gene.fa --probes probes.txt"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import core


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="transquant-w",
        description=(
            "Compute the TransQuant probe weight factor W, gene length L and probe "
            "localisation profile N for an smFISH probe library "
            "(Bahar Halpern & Itzkovitz, Methods 2016, Eq. 7)."
        ),
    )
    p.add_argument("--target", required=True, type=Path,
                   help="FASTA (or plain text) of the FULL genomic locus, in the direction of transcription.")
    p.add_argument("--probes", required=True, type=Path,
                   help="Probe sequences: one per line, or a tab-delimited table with a sequence column.")
    p.add_argument("--out", type=Path, default=None,
                   help="Output prefix. Writes <out>.png; with --profile-tsv also <out>_profile.tsv. "
                        "Default: no files, results printed only.")
    p.add_argument("--svg", action="store_true", help="Also write <out>.svg.")
    p.add_argument("--profile-tsv", action="store_true",
                   help="Write N(i) for every position to <out>_profile.tsv.")
    p.add_argument("--no-revcomp", action="store_true",
                   help="Probes are already target-sense; do not reverse complement them.")
    p.add_argument("--no-concat", action="store_true",
                   help="For multi-entry FASTA, use only the first entry instead of concatenating.")
    p.add_argument("--no-case-shading", action="store_true",
                   help="Do not shade UPPERCASE runs of the target (UCSC exons) on the RNA bar of the plot.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        target_text = args.target.read_text()
        target, w1 = core.parse_target(target_text, concatenate=not args.no_concat)
        probes, w2 = core.parse_probes(args.probes.read_text())
        res = core.compute(target, probes, revcomp_probes=not args.no_revcomp, warnings=w1 + w2)
    except (OSError, core.InputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    for w in res.warnings:
        print(f"warning: {w}", file=sys.stderr)

    print(f"W\t{res.W:.4f}")
    print(f"L\t{res.L}")
    print(f"probes_supplied\t{res.n_probes}")
    print(f"probes_matched\t{res.n_matched_probes}")
    print(f"binding_sites\t{res.n_sites}")

    if args.out is not None:
        blocks = None if args.no_case_shading else core.case_blocks(target_text, concatenate=not args.no_concat)
        if blocks:
            print("note: uppercase runs of the target are shaded on the RNA bar (UCSC writes exons in "
                  "uppercase; check this matches your download). Disable with --no-case-shading.", file=sys.stderr)
        fig = core.make_plot(res, blocks=blocks)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        Path(f"{args.out}.png").write_bytes(core.figure_bytes(fig, "png"))
        if args.svg:
            Path(f"{args.out}.svg").write_bytes(core.figure_bytes(fig, "svg"))
        if args.profile_tsv:
            with open(f"{args.out}_profile.tsv", "w") as fh:
                fh.write("position\tN\n")
                fh.writelines(f"{i}\t{n}\n" for i, n in enumerate(res.N, start=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
