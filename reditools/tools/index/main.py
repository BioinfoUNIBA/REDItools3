"""Commandline tool for REDItools."""

import sys
from reditools.file_utils import open_stream
from reditools.region import Region
from reditools.rtindexer import RTIndexer
from reditools.tools.index.parse_args import parse_args


def main():
    """Perform RNA editing analysis."""
    options = parse_args()
    if options.region:
        indexer = RTIndexer(Region.parse_string(options.region))
    else:
        indexer = RTIndexer()

    if options.exclude_regions:
        for _ in options.exclude_regions:
            indexer.add_exclusions_from_bed(_)

    if options.bed_file:
        for _ in options.bed_file:
            indexer.add_target_from_bed(_)

    if options.output_file:
        stream = open_stream(options.output_file, 'w')
    else:
        stream = sys.stdout

    for _ in options.file:
        indexer.add_rt_output(_)

    for nuc, idx in sorted(indexer.calc_index().items()):
        stream.write(f'{nuc}\t{idx}\n')


if __name__ == '__main__':
    main()
