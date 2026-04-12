import csv
from tempfile import NamedTemporaryFile
from reditools.logger import Logger

_empty = '-'

def mean_quality(rtresult):
    if len(rtresult) == 0:
        return 0
    return sum(rtresult.qualities) / len(rtresult)

def edit_ratio(rtresult):
    max_edits = 0
    for base, count in zip(rtresult._bases, rtresult):
        if base != rtresult.reference and count > max_edits:
            max_edits = count
    try:
        return max_edits / (rtresult['REF'] + max_edits)
    except ZeroDivisionError:
        return 0

def write_results(rtresults, file_name, output_format,
                  temp_dir, filters, logger):
    """
    Write the results from a REDItools analysis to a temporary file.

    Parameters:
        rtresults (iterable): REDItools results
        file_name (string): Input file name for analysis
        output_format (dict): keyword arguments for csv.writer constructor.
        temp_dir (str): Location to save results
        filters (RTChecks): Filters what to write

    Returns:
        string: Name of the results file.
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
                f'{mean_quality(rt_result):.2f}',
                list(rt_result),
                ' '.join(sorted(variants)) if variants else _empty,
                f'{edit_ratio(rt_result):.2f}',
                _empty, _empty, _empty, _empty, _empty,
            ])
        return stream.name
