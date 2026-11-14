"""CLI arguments for legacy options."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools import reditools
from reditools.tools.analyze.parse_args.bounded_types import bounded_int

if TYPE_CHECKING:
    import argparse

def add_leg_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for legacy options.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    leg_group = parser.add_argument_group(
        title="Legacy Options",
    )
    leg_group.add_argument(
        "-e",
        "--exclude-multis",
        default=False,
        help=(
            "Do not report any position with more than one alternate base. "
            "(Equivalent to -Men/--max-editing-nucleotides 1)"
        ),
        action="store_true",
    )
    leg_group.add_argument(
        "-N",
        "--dna",
        default=False,
        help=(
            "Run REDItools on DNA-Seq data. "
            f"(Equivalent to -s/--strand {reditools.UNSTRANDED_MODE})"
        ),
        action="store_true",
    )
    leg_group.add_argument(
        "-m",
        "--load-omopolymeric-file",
        help=(
            "BED file of homopolymeric positions. Regions in the BED file "
            "will be excluded from the analysis. (Same effect as providing "
            "the BED file with the -k/--exclude-regions option."
        ),
    )
    leg_group.add_argument(
        "-S",
        "--strict",
        help=(
            "Activate strict mode: only sites with edits will be included in "
            "the output. (Equivalent to -me/--min-edits 1)"
        ),
        action="store_true",
    )
    leg_group.add_argument(
        "-sf",
        "--splicing-file",
        help=(
            "The splicing file is a space delimited file with five columns: "
            "chromosome, start position (one-index inclusive), stop "
            "(ignored),  splice (either A for acceptor or D for donor), and "
            "strand (either + or -). A header is optional, but must start "
            "with #. Used in conjunctions with -ss/--splicing-span."
        ),
    )
    leg_group.add_argument(
        "-ss",
        "--splicing-span",
        type=bounded_int(min_value=1),
        default=4,
        help=(
            "The splicing span. Used in conjunction with -sf/--splicing-file."
        ),
    )
