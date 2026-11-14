"""CLI arguments for filtering results."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools.tools.analyze.parse_args.bounded_types import (
    bounded_float,
    bounded_int,
)

if TYPE_CHECKING:
    import argparse

def add_rf_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for filtering results.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    rf_group = parser.add_argument_group(
        title="Result Filters",
    )
    rf_group.add_argument(
        "-men",
        "--min-edits-per-nucleotide",
        type=int,
        default=0,
        help=(
            "Position where any variant has a frequency less than -men but "
            "more than zero will not be reported. (Corresponds to the "
            "BaseCount column.)"
        ),
    )
    rf_group.add_argument(
        "-me",
        "--min-edits",
        type=bounded_int(0, 4),
        default=1,
        help=(
            "Positions with fewer than -me unique variants (listed in the "
            "AllSubs column) will be excluded from the results."
        ),
    )
    rf_group.add_argument(
        "-Men",
        "--max-editing-nucleotides",
        type=bounded_int(min_value=0, max_value=4),
        default=4,
        help=(
            "Positions with more than -Men unique variants (listed in the "
            "AllSubs column) will be excluded from the results."
        ),
    )
    rf_group.add_argument(
        "-v",
        "--variants",
        nargs="*",
        default=["all"],
        help=(
            "Which editing events to report. Each edit should be two "
            "characters and separated by spaces (e.g. AG CT). Use 'all' to "
            "report all variants. (Corresponds to the AllSubs column)"
        ),
    )
    rf_group.add_argument(
        "-l",
        "--min-read-depth",
        type=bounded_int(min_value=1),
        default=1,
        help=(
            "Only report on positions with at least -l reads (corresponds to "
            "the Coverage column.)"
        ),
    )
    rf_group.add_argument(
        "-n",
        "--noise-threshold",
        type=bounded_float(min_value=0, max_value=1),
        default=0,
        help=(
            "Only report edits with a ratio greater than or equal to -n "
            "(corresponds to the Frequency column.) This option generally "
            "filters what gets reported in the AllSubs column. Edits filtered "
            "by this option will not be considered by --variants, "
            "--max-editing-nucleotides, --min-edits, or "
            "--min-edits-per-nucleotide filters."
        ),
    )
