"""CLI options relating to base quality filters."""
from __future__ import annotations

import argparse

from reditools.tools.analyze.parse_args.bounded_types import bounded_int


def add_bqf_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for base quality filters.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    bqf_group = parser.add_argument_group(
        title="Base/Read Quality Controls",
    )
    bqf_group.add_argument(
        "-mrl",
        "--min-read-length",
        type=int,
        default=30,
        help="reads shorter than -mrl will be discarded.",
    )
    bqf_group.add_argument(
        "-q",
        "--min-read-quality",
        type=int,
        default=20,
        help="Reads with MAPQ below -q will be discarded.",
    )
    bqf_group.add_argument(
        "-bq",
        "--min-base-quality",
        type=int,
        default=30,
        help="Bases with a Phred quality score below -bq will bed discarded.",
    )
    bqf_group.add_argument(
        "-mbp",
        "--min-base-position",
        type=bounded_int(min_value=0),
        default=0,
        help="Ignores the first -mbp bases in each read.",
    )
    bqf_group.add_argument(
        "-Mbp",
        "--max-base-position",
        type=bounded_int(min_value=0),
        default=0,
        help="Ignores the last -Mpb bases in each read.",
    )
    bqf_group.add_argument(
        "-E",
        "--exclude-reads",
        help="Text file listing read names to exclude from analysis.",
    )
    bqf_group.add_argument(
        "--exclude_reads",
        help=argparse.SUPPRESS,
    )
