"""Commandline tool for REDItools."""

import argparse
import csv
import sys
from itertools import permutations
from json import loads as load_json

from reditools.file_utils import open_stream, read_bed_file
from reditools.region import Region

_ref = 'Reference'
_position = 'Position'
_contig = 'Region'
_count = 'BaseCount[A,C,G,T]'
_strand = 'Strand'
_nucs = 'ACGT'
_ref_set = {f'{nuc}-{nuc}' for nuc in _nucs}


class Index(object):
    """Utility for calculating editing indices."""

    def __init__(self, region=None, strand=0):
        """
        Create a new Index.

        Parameters:
            region (Region): Limit results to the given genomic region
            strand (int): Either 0, 1, or 2 for unstranded, reverse, or forward
        """
        self.targets = {}
        self.exclusions = {}
        self.counts = {'-'.join(_): 0 for _ in permutations(_nucs, 2)}
        self.region = region
        self.strand = ['*', '-', '+'][strand]

    def add_target_from_bed(self, fname):
        """
        Only report index data for regions from a given bed file.

        Parameters:
            fname (str): Path to BED formatted file.
        """
        for region in read_bed_file(fname):
            self.targets[region.contig] = update_region_dict(
                self.targets,
                region,
            )

    def add_exclusions_from_bed(self, fname):
        """
        Exclude index data for regions from a given bed file.

        Parameters:
            fname (str): Path to BED formatted file.
        """
        for region in read_bed_file(fname):
            self.exclusions[region.contig] = update_region_dict(
                self.exclusions,
                region,
            )

    def in_region_list(self, region_list, contig, position):
        """
        Check if a genomic position is in a list of regions.

        Parameters:
            region_list (dict): Region list to check
            contig (str): Contig/Chromsome name
            position (int): Coordinate

        Returns:
            True if the position is present, else False
        """
        return position in region_list.get(contig, [])

    def in_targets(self, contig, position):
        """
        Check if a genomic position is in the target list.

        Parameters:
            contig (str): Contig/Chromsome name
            position (int): Coordiante

        Returns:
            True if there are no targets or the position is in the target
            list; else False
        """
        return not self.targets or self.in_region_list(self.targets)

    def in_exclusions(self, contig, position):
        """
        Check if a genomic position is in the exclusions list.

        Parameters:
            contig (str): Contig/Chromsome name
            position (int): Coordiante

        Returns:
            True if there are no exclusions or the position is in the
            exclusions list; else False
        """
        return self.exclusions and self.in_region_list(self.exclusions)

    def do_ignore(self, row):
        """
        Check whether a row should meets analysis criteria.

        Parameters:
            row (dict): Row from REIDtools output file.

        Returns:
            True if the row should be discarded; else False
        """
        if '*' != self.strand != row[_strand]:
            return True
        if self.region:
            if not self.region.contains(row[_contig], row[_position]):
                return True
        if self.in_exclusions(row[_contig], row[_position]):
            return True
        return not self.in_targets(row[_contig], row[_position])

    def add_rt_output(self, fname):
        """
        Count the number of reads with matches and substitutions.

        Parameters:
            fname (str): File path to a REDItools output
        """
        stream = open_stream(fname)
        reader = csv.DictReader(stream, delimiter='\t')
        for row in reader:
            if self.do_ignore(row):
                continue
            ref = row[_ref]
            reads = load_json(row[_count])
            for nuc, count in zip(_nucs, reads):
                key = f'{nuc}-{ref}'
                self.counts[key] = self.counts.get(key, 0) + count
        stream.close()

    def calc_index(self):
        """
        Compute all editing indices.

        Returns:
            Dictionary of indices
        """
        keys = set(self.counts) - _ref_set
        indices = {}
        for idx in keys:
            ref = idx[-1]
            numerator = self.counts[idx]
            denominator = self.counts.get(self.ref_edit(ref), 0) + numerator
            if denominator == 0:
                indices[idx] = 0
            else:
                indices[idx] = numerator / denominator
        return indices

    def ref_edit(self, ref):
        """
        Format a base as a non-edit.

        Parameters:
            ref (str): Reference base

        Returns:
            A string in the format of {ref}-{ref}
        """
        return f'{ref}-{ref}'


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


def main():
    """Perform RNA editing analysis."""
    options = parse_options()
    if options.region:
        indexer = Index(Region(string=options.region), strand=options.strand)
    else:
        indexer = Index(strand=options.strand)

    if options.exclude_regions:
        for exc_fname in options.exclude_regions:
            indexer.add_exclusions_from_bed(exc_fname)

    if options.bed_file:
        for trg_fname in options.bed_file:
            indexer.add_target_from_bed(trg_fname)

    if options.output_file:
        stream = open_stream(options.output_fipe, 'w')
    else:
        stream = sys.stdout

    for fname in options.file:
        indexer.add_rt_output(fname)

    for nuc, idx in sorted(indexer.calc_index().items()):
        stream.write(f'{nuc}\t{idx}\n')


def update_region_dict(region_dict, region):
    """
    Add a region to a region dictionary.

    Parameters:
        region_dict (dict): Region dictionary
        region (Region): Region to add

    Returns:
        An updated copy of region_dict
    """
    return region_dict.get(region.contig, set()) | region.enumerate()


if __name__ == '__main__':
    main()
