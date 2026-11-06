"""Write analysis results."""

import csv
from pathlib import Path
from typing import Callable, Iterable

from reditools.compiled_position import RTResult
from reditools.logger import Logger
from reditools.tools.analyze.rtchecks import RTChecks


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
    with Path(filename).open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        for rt_result in rtresults:
            msg = filters.check(rt_result)
            if msg:
                logger(Logger.debug_level, *msg)
                continue
            writer.writerow(rt_result.output_list(strand_numbers))
