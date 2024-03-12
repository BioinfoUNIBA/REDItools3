"""Repeat Sequence Identifier."""
import argparse
import sys

from pysam import FastaFile

from reditools import file_utils


def find_homo_seqs(seq, length=5):
    """
    Locate regions of repeated bases.

    Parameters:
        seq (str): The DNA sequence
        length (int): Minimum number of sequential repeats.

    Yields:
        start, stop, base
    """
    h_base = None
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


def parse_options():
    """
    Parse commandline arguments.

    Returns:
        namespace
    """
    parser = argparse.ArgumentParser(description='REDItools 2.0')
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
        help='Destination to write results. Default is to use STDOUT. ' +
        'If the filename ends in .gz, the contents will be gzipped.',
    )

    return parser.parse_args()


def main():
    """Report repetative regions."""
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

    for seq_name in fasta.references:
        seq = fasta.fetch(seq_name)
        for region in find_homo_seqs(seq, options.min_length):
            fields = [
                seq_name,
                region[0],
                region[1],
                region[1] - region[0],
                region[2],
            ]
            as_str = [str(_) for _ in fields]
            stream.write('\t'.join(as_str) + '\n')
