
import csv
import sys

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
    tfs: list[str],
    output_file: str | None = None,
    mode: str = 'w',
    encoding: str = 'utf-8',
    **format_args,
) -> None:
    """Concatenate temporary results files into the final output.

    Parameters
    ----------
    tfs : list[str]
        List of paths to temporary files.
    output_file : str, optional
        Path to the final output file. If None, write to stdout.
    mode : str, default 'w'
        Mode in which the output file is opened.
    encoding : str, default 'utf-8'
        Encoding for the output file.
    **format_args
        Additional arguments for csv.writer.
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
