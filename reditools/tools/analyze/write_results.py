"""Write analysis results."""

import csv
from pathlib import Path
from typing import Callable, Iterable

from reditools.compiled_position import RTResult
from reditools.constants import (
    forward_strand_number,
    forward_strand_symbol,
    reverse_strand_number,
    reverse_strand_symbol,
    undetermined_strand_number,
    undetermined_strand_symbol,
)
from reditools.logger import Logger
from reditools.tools.analyze.rtchecks import RTChecks

_empty = "-"

def write_results(
        rtresults: Iterable[RTResult],
        filename: str,
        filters: RTChecks,
        logger: Callable,
        strand_numbers: bool = False,
) -> None:
    """Write analysis results to a file.

    Parameters
    ----------
    rtresults : Iterable[RTResult]
        The analysis results for each position.
    filename : str
        Where to save the results.
    filters : RTChecks
        The quality control checks to apply.
    logger : Callable
        The logger function for debug messages.
    strand_numbers : bool
        If False, uses -, +, and * for the Strand column.
        If True, uses 0, 1, and 2 for the Strand column.
    """
    if strand_numbers:
        strand_lookup = {
            forward_strand_symbol: forward_strand_number,
            reverse_strand_symbol: reverse_strand_number,
            undetermined_strand_symbol: undetermined_strand_number,
        }
    else:
        strand_lookup = {
            forward_strand_symbol: forward_strand_symbol,
            reverse_strand_symbol: reverse_strand_symbol,
            undetermined_strand_symbol: undetermined_strand_symbol,
        }
    with Path(filename).open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        for rt_result in rtresults:
            msg = filters.check(rt_result)
            if msg:
                logger(Logger.debug_level, *msg)
                continue
            writer.writerow([
                rt_result.contig,
                rt_result.position + 1,
                rt_result.reference,
                strand_lookup[rt_result.strand],
                len(rt_result),
                f"{rt_result.mean_quality:.2f}",
                list(rt_result),
                " ".join(rt_result.variants) if rt_result.variants else _empty,
                f"{rt_result.edit_ratio:.2f}",
                _empty, _empty, _empty, _empty, _empty,
            ])
