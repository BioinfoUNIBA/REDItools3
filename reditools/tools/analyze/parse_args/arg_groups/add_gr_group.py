"""CLI arguments for filtering by genomic region."""
from __future__ import annotations

import argparse


def add_gr_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for filtering by genomic region.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    gr_group = parser.add_argument_group(
        title="Genomic Region Filters",
    )
    gr_group.add_argument(
        "-k",
        "--exclude-regions",
        nargs="+",
        help="Do not report on regions in the provided BED file(s).",
    )
    gr_group.add_argument(
        "--exclude_regions",
        nargs="+",
        help=argparse.SUPPRESS,
    )
    gr_group.add_argument(
        "-B",
        "--bed-file",
        nargs="+",
        help="Only reports on regions in the provided BED file(s).",
    )
    gr_group.add_argument(
        "--bed_file",
        nargs="+",
        help=argparse.SUPPRESS,
    )
