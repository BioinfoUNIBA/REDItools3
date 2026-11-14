"""CLI arguments for parallel processing."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools.tools.analyze.parse_args.bounded_types import bounded_int

if TYPE_CHECKING:
    import argparse

def add_para_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for parallel processing.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    para_group = parser.add_argument_group(
        title="Parallel Processing Options",
    )
    para_group.add_argument(
        "-t",
        "--threads",
        help=(
            "Number of threads for parallel processing. Note that the "
            "maximum number of usable threads is equivalent to the number of "
            "chromosomes in your alignment genome unless you use the --window "
            "option."
        ),
        type=bounded_int(min_value=1),
        default=1,
    )
    para_group.add_argument(
        "-w",
        "--window",
        help=(
            "How many bp should be processed by each thread at a time. "
            "Zero uses the full contig."
        ),
        type=bounded_int(min_value=0),
        default=0,
    )
