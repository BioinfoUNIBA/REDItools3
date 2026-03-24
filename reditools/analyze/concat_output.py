"""Commandline tool for REDItools."""

import sys
import csv
from reditools import file_utils

fieldnames = [
    'Region',
    'Position',
    'Reference',
    'Strand',
    'Coverage',
    'MeanQ',
    'BaseCount[A,C,G,T]',
    'AllSubs',
    'Frequency',
    'gCoverage',
    'gMeanQ',
    'gBaseCount[A,C,G,T]',
    'gAllSubs',
    'gFrequency',
]


def concat_output(
    tfs,
    output_file=None,
    mode='w',
    encoding='utf-8',
    **format_args,
):
    """
    Write the output of a REDItools analysis.

    Parameters:
        options (namespace): Commandline options for file formatting.
        tfs (list): Temporary files containing REDItools results
    """
    # Setup final output file
    if output_file is None:
        stream = sys.stdout
    else:
        stream = file_utils.open_stream(
            output_file,
            mode,
            encoding,
        )

    with stream:
        writer = csv.writer(stream, **format_args)
        if 'a' not in mode:
            writer.writerow(fieldnames)
        file_utils.concat(stream, *tfs, encoding=encoding)
