"""Commandline tool for REDItools."""

import argparse
import csv
from json import loads as load_json

_ref = 'Reference'
_count = 'BaseCount[A,C,G,T]'
_nucs = 'ACGT'
_ref_set = {f'{nuc}-{nuc}' for nuc in _nucs}


def count_nuc_events(fname):
    """
    Count the number of reads with matches and substitutions.

    Parameters:
        fname (str): File path to a REDItools output

    Returns:
        Dict with keys in the format of base-base
    """
    counts = {}
    with open(fname, 'r') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        for row in reader:
            ref = row[_ref]
            reads = load_json(row[_count])
            for nuc, count in zip(_nucs, reads):
                key = f'{nuc}-{ref}'
                counts[key] = counts.get(key, 0) + count
    return counts


def merge(dict_list):
    """
    Combine multiple dictionaries together using the sum of their values.

    Parameters:
        dict_list (list): A list of dictionariees to merge

    Returns:
        A final cumulative dictionary
    """
    if len(dict_list) == 1:
        return dict_list[0]
    merged = merge(dict_list[1:])
    for key, count in dict_list[0].items():
        merged[key] = merged.get(key, 0) + count
    return merged


def index(rt_output_files):
    """
    Compute all editing indices.

    Parameters:
        rt_output_files (list): All REDItools output files to be processed
    """
    counts = merge([count_nuc_events(fname) for fname in rt_output_files])
    indices = set(counts) - _ref_set
    for idx in indices:
        ref = idx[-1]
        numerator = counts[idx]
        denominator = counts[f'{ref}-{ref}'] + numerator
        scale_factor = numerator / denominator
        print(f'{idx}\t{scale_factor}')


def parse_options():  # noqa:WPS213
    """
    Parse commandline options for REDItools.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument(
        'file',
        nargs='+',
        help='The REDItools output file to be analyzed',
    )
    parser.add_argument(
        '-o',
        '--output-file',
        help='The output statistics file',
    )
    parser.add_argument(
        '-s',
        '--strand',
        choices=(0, 1, 2),
        type=int,
        default=0,
        help='Strand: this can be 0 (unstranded),' +
        '1 (secondstrand oriented) or ' +
        '2 (firststrand oriented)',
    )
    parser.add_argument(
        '-g',
        '--region',
        help='The region of the bam file to be analyzed',
    )
    parser.add_argument(
        '-B',
        '--bed_file',
        help='Path of BED file containing target regions',
    )
    parser.add_argument(
        '-k',
        '--exclude_regions',
        nargs='+',
        help='Path of BED file containing regions to exclude from analysis',
    )

    return parser.parse_args()


def main():
    """Perform RNA editing analysis."""
    options = parse_options()
    index(options.file[0])


if __name__ == '__main__':
    main()
