import argparse
import csv
import sys
from typing import Iterator

from pysam import FastaFile

from reditools import file_utils


def find_homo_seqs(seq: str, length: int = 5) -> Iterator[tuple[int, int, str]]:
    """Find homopolymeric sequences in a string.

    Parameters
    ----------
    seq : str
        The sequence to search.
    length : int, optional
        The minimum length of the repeat, by default 5.

    Yields
    ------
    Iterator[tuple[int, int, str]]
        A tuple containing (start, end, base) for each homopolymeric sequence.
    """
    h_base = ''
    start = 0
    count = 0

    for pos, base in enumerate(seq):
        if base == h_base:
            count += 1
        else:
            if count >= length:
                yield (start, start + count, h_base)
            count = 0
            start = pos
            h_base = base
    if count >= length:
        yield (start, start + count, h_base)


def parse_options() -> argparse.Namespace:
    """Parse command-line arguments for reditools find-repeats.

    Returns
    -------
    argparse.Namespace
        The parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="reditools find-repeats",
        description='REDItools3',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'file',
        help='The fasta file to be analyzed',
    )
    parser.add_argument(
        '-l',
        '--min-length',
        type=int,
        default=5,
        help='Minimum length of repeat region',
    )
    parser.add_argument(
        '-o',
        '--output',
        default='/dev/stdout',
        help='Destination to write results. Default is to use STDOUT. ' +
        'If the filename ends in .gz, the contents will be gzipped.',
    )

    return parser.parse_args()

def iter_homo_output(
        fasta: FastaFile,
        min_length: int,
) -> Iterator[tuple[str, int, int, int, str]]:
    """Iterate over homopolymeric sequences in a FASTA file.

    Parameters
    ----------
    fasta : FastaFile
        The opened FASTA file.
    min_length : int
        The minimum length of the repeat.

    Yields
    ------
    Iterator[tuple[str, int, int, int, str]]
        A tuple containing (seq_name, start, end, length, base) for each repeat.
    """
    for seq_name in fasta.references:
        for region in find_homo_seqs(fasta.fetch(seq_name), min_length):
            yield (
                seq_name,
                region[0],
                region[1],
                region[1] - region[0],
                region[2],
            )

def main() -> None:
    """Execute the reditools find-repeats tool.

    This tool identifies homopolymeric sequences in a FASTA file.
    """
    options = parse_options()
    fasta = FastaFile(options.file)

    if options.output:
        stream = file_utils.open_stream(
            options.output,
            'wt',
            encoding='utf-8',
        )
    else:
        stream = sys.stdout

    writer = csv.writer(stream, delimiter='\t')
    writer.writerows(iter_homo_output(fasta, options.min_length))
