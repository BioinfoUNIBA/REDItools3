import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for reditools index.

    Returns
    -------
    argparse.Namespace
        The parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="reditools index",
        description='REDItools3',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'file',
        nargs='+',
        help='The REDItools output file to be analyzed',
    )
    parser.add_argument(
        '-o',
        '--output-file',
        default='/dev/stdout',
        help='The output statistics file',
    )
    parser.add_argument(
        '-g',
        '--region',
        help='The genomic region to be analyzed',
    )
    parser.add_argument(
        '-B',
        '--bed_file',
        nargs='+',
        help='Path of BED file containing target regions',
    )
    parser.add_argument(
        '-k',
        '--exclude_regions',
        nargs='+',
        help='Path of BED file containing regions to exclude from analysis',
    )

    return parser.parse_args()
