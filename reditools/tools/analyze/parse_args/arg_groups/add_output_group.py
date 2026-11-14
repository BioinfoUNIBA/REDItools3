"""CLI arguments for output location and formatting."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse

def add_output_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for output location and formatting.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    output_group = parser.add_argument_group(
        title="Output Options",
    )
    output_group.add_argument(
        "-o",
        "--output-file",
        help="Path to write output to.",
        default="/dev/stdout",
    )
    output_group.add_argument(
        "-a",
        "--append-file",
        action="store_true",
        help="Appends results to file (and creates if not existing).",
    )
    output_group.add_argument(
        "--number-strand-output",
        action="store_true",
        help=(
            "Replaces the '-', '+', and '*' values with 0, 1, and 2 in the "
            "Strand column of the output."
        ),
    )
