import csv
from tempfile import NamedTemporaryFile
from typing import Callable, Iterator

from reditools.compiled_position import RTResult
from reditools.logger import Logger
from reditools.tools.analyze.rtchecks import RTChecks

_empty = '-'

def write_results(
        rtresults: Iterator[RTResult],
        output_format: dict,
        temp_dir: str,
        filters: RTChecks,
        logger: Callable,
) -> str:
    """Write analysis results to a temporary file.

    Parameters
    ----------
    rtresults : Iterator[RTResult]
        The analysis results for each position.
    output_format : dict
        The CSV formatting options.
    temp_dir : str
        The directory to save temporary files.
    filters : RTChecks
        The quality control checks to apply.
    logger : Callable
        The logger function for debug messages.

    Returns
    -------
    str
        The path to the temporary file containing the results.
    """
    with NamedTemporaryFile(mode='w', delete=False, dir=temp_dir) as stream:
        writer = csv.writer(stream, **output_format)
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
                f'{rt_result.mean_quality:.2f}',
                list(rt_result),
                ' '.join(sorted(variants)) if variants else _empty,
                f'{rt_result.edit_ratio:.2f}',
                _empty, _empty, _empty, _empty, _empty,
            ])
        return stream.name
