"""Write analysis results."""

import csv
from pathlib import Path
from typing import Callable, Iterator

from reditools.compiled_position import RTResult
from reditools.logger import Logger
from reditools.tools.analyze.rtchecks import RTChecks

_empty = "-"

def write_results(
        rtresults: Iterator[RTResult],
        filename: str,
        filters: RTChecks,
        logger: Callable,
) -> None:
    """Write analysis results to a file.

    Parameters
    ----------
    rtresults : Iterator[RTResult]
        The analysis results for each position.
    filename : str
        Where to save the results.
    filters : RTChecks
        The quality control checks to apply.
    logger : Callable
        The logger function for debug messages.
    """
    with Path(filename).open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        for rt_result in rtresults:
            msg = filters.check(rt_result)
            if msg:
                logger(Logger.debug_level, *msg)
                continue
            variants = rt_result.variants
            writer.writerow([
                rt_result.contig,
                rt_result.position + 1,
                rt_result.reference,
                rt_result.strand,
                len(rt_result),
                f"{rt_result.mean_quality:.2f}",
                list(rt_result),
                " ".join(sorted(variants)) if variants else _empty,
                f"{rt_result.edit_ratio:.2f}",
                _empty, _empty, _empty, _empty, _empty,
            ])
