"""CLI arguments for DNA-RNA combo analysis."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools.tools.analyze.parse_args.bounded_types import bounded_int

if TYPE_CHECKING:
    import argparse

def add_dna_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for DNA-RNA combo analysis.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    dna_group = parser.add_argument_group(
        title="Annotate Options",
        description=(
            "Use these options to analyze a DNA sample in parallel, using it "
            "to fill in the last five columns of the output. "
            "This essentially combines the analyze tool with the annotate "
            "tool."
        ),
    )
    dna_group.add_argument(
        "-D",
        "--dna-sample",
        help="BAM file containing DNA alignment.",
    )
    dna_group.add_argument(
        "-Dmrl",
        "--dna-min-read-length",
        type=int,
        default=30,
        help="DNA reads shorter than -mrl will be discarded.",
    )
    dna_group.add_argument(
        "-Dq",
        "--dna-min-read-quality",
        type=int,
        default=20,
        help="DNA reads with MAPQ below -q will be discarded.",
    )
    dna_group.add_argument(
        "-Dbq",
        "--dna-min-base-quality",
        type=int,
        default=30,
        help=(
            "DNA bases with a Phred quality score below -bq will bed discarded."
        ),
    )
    dna_group.add_argument(
        "-Dmbp",
        "--dna-min-base-position",
        type=bounded_int(min_value=0),
        default=0,
        help="Ignores the first -mbp bases in each DNA read.",
    )
    dna_group.add_argument(
        "-DMbp",
        "--dna-max-base-position",
        type=bounded_int(min_value=0),
        default=0,
        help="Ignores the last -Mpb bases in each DNA read.",
    )
    dna_group.add_argument(
        "-Dl",
        "--dna-min-read-depth",
        type=bounded_int(min_value=0),
        default=0,
        help=(
            "Only annotate positions with at least -Dl DNA reads (corresponds "
            "the gCoverage column.)"
        ),
    )
    dna_group.add_argument(
        "-Df",
        "--dna-filter",
        help="Only report on positions with DNA coverage.",
        action="store_true",
    )
