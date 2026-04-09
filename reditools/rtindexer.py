import csv
from itertools import permutations
from json import loads as load_json

from reditools.file_utils import open_stream, read_bed_file
from reditools.region_collection import RegionCollection

_ref = 'Reference'
_position = 'Position'
_contig = 'Region'
_count = 'BaseCount[A,C,G,T]'
_strand = 'Strand'
_nucs = 'ACGT'


class RTIndexer(object):
    """Utility for calculating editing indices."""

    def __init__(self, region=None, strand=0):
        """
        Create a new Index.

        Parameters:
            region (Region): Limit results to the given genomic region
            strand (int): Either 0, 1, or 2 for unstranded, reverse, or forward
        """
        self.targets = False
        self.exclusions = False
        self.counts = {'-'.join(_): 0 for _ in permutations(_nucs, 2)}
        self.region = region
        self.strand = ['*', '-', '+'][strand]

    def add_target_from_bed(self, fname):
        """
        Only report index data for regions from a given bed file.

        Parameters:
            fname (str): Path to BED formatted file.
        """
        if not self.targets:
            self.targets = RegionCollection()
        self.targets.add_regions(read_bed_file(fname))

    def add_exclusions_from_bed(self, fname):
        """
        Exclude index data for regions from a given bed file.

        Parameters:
            fname (str): Path to BED formatted file.
        """
        if not self.exclusions:
            self.exclusions = RegionCollection()
        self.exclusions.add_regions(read_bed_file(fname))

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
        return not self.targets or self.targets.contains(contig, position)

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
        return self.exclusions and self.exclusions.contains(contig, position)

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
            position = int(row[_position])
            if self.region[0] != row[_contig] or \
                    self.region[1] > position or \
                    self.region[2] is not None and self.region[2] < position:
                return True
        if self.in_exclusions(row[_contig], int(row[_position])):
            return True
        return not self.in_targets(row[_contig], int(row[_position]))

    def add_rt_output(self, fname):
        """
        Count the number of reads with matches and substitutions.

        Parameters:
            fname (str): File path to a REDItools output
        """
        with open_stream(fname) as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if self.do_ignore(row):
                    continue
                for nuc, count in zip(_nucs, load_json(row[_count])):
                    key = f'{nuc}-{row[_ref]}'
                    self.counts[key] = self.counts.get(key, 0) + count

    def calc_index(self):
        """
        Compute all editing indices.

        Returns:
            Dictionary of indices
        """
        indices = {}
        for idx in set(self.counts) - {f'{nuc}-{nuc}' for nuc in _nucs}:
            ref = idx[-1]
            numerator = self.counts[idx]
            denominator = self.counts.get(self.ref_edit(ref), 0) + numerator
            if denominator == 0:
                indices[idx] = 0
            else:
                indices[idx] = 100 * numerator / denominator
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
