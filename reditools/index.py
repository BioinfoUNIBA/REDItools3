"""Commandline tool for REDItools."""

import argparse
import csv

_ref = 'Reference'
_count = 'BaseCount[A,C,G,T]'
_bases = 'ACGT'

def go(rt_output_file):
    counts = {}
    with open(rt_output_file, "r") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            ref = row[_ref]
            reads = json.loads(row[_count])
            for base, count in zip(_bases, reads):
                key=f"{base}-{ref}"
                counts[key] = counts.get(key, 0) + count
    indices = set(counts) - set([f"{base}-{base}" for base in _bases])
    for idx in indices:
        ref = idx[-1]
        numerator = counts[idx]
        denominator = counts[f"{ref}-{ref}"] + numerator
        print(f"{idx}\t{numerator / denominator}")

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
        help='The self.region of the bam file to be analyzed',
    )
   parser.add_argument(
        '-B',
        '--bed_file',
        help='Path of BED file containing target self.regions',
    )
   parser.add_argument(
        '-k',
        '--exclude_regions',
        nargs='+',
        help='Path of BED file containing regions to exclude from analysis',
    )
    parser.add_argument(
        '-d',
        '--debug',
        default=False,
        help='REDItools is run in DEBUG mode.',
        action='store_true',
    )

    return parser.parse_args()


def main():
    """Perform RNA editing analysis."""
    options = parse_options()
    go(options.file)

if __name__ == "__main__":
    main()
