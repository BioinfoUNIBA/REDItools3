import csv
from itertools import permutations

from reditools.file_utils import open_stream, read_bed_file
from reditools.region_collection import RegionCollection

_ref = 'Reference'
_position = 'Position'
_contig = 'Region'
_count = 'BaseCount[A,C,G,T]'
_strand = 'Strand'
_nucs = 'ACGT'


class RTIndexer(object):
    """
    Calculate editing indices from REDItools output.

    Parameters
    ----------
    region : tuple[str, int, int | None] | None, optional
        Genomic region (contig, start, stop) to limit analysis (default is
        None).
    """

    def __init__(
            self,
            region: tuple[str, int, int | None] | None=None,
    ):
        """
        Initialize the RTIndexer.

        Parameters
        ----------
        region : tuple[str, int, int | None] | None, optional
            Genomic region (contig, start, stop) to limit analysis (default is
            None).
        """
        self.targets = RegionCollection()
        self.exclusions = RegionCollection()
        self.counts = {'-'.join(_): 0 for _ in permutations(_nucs, 2)}
        self.region = region

    def add_target_from_bed(self, fname: str) -> None:
        """
        Add target regions from a BED file.

        Parameters
        ----------
        fname : str
            Path to the BED file.
        """
        self.targets.add_regions(read_bed_file(fname))

    def add_exclusions_from_bed(self, fname: str) -> None:
        """
        Exclude regions from a BED file.

        Parameters
        ----------
        fname : str
            Path to the BED file.
        """
        self.exclusions.add_regions(read_bed_file(fname))

    def do_ignore(self, row: dict) -> bool:
        """
        Check if a row from REDItools output should be ignored.

        Parameters
        ----------
        row : dict
            A dictionary representing a row of REDItools output.

        Returns
        -------
        bool
            True if the row should be ignored, False otherwise.
        """
        if self.region:
            position = int(row[_position])
            if self.region[0] != row[_contig] or \
                    self.region[1] > position or \
                    self.region[2] is not None and self.region[2] < position:
                return True
        if self.exclusions and self.exclusions.contains(
                row[_contig],
                int(row[_position]),
        ):
            return True
        if self.targets:
             return not self.targets.contains(
                 row[_contig],
                 int(row[_position]),
            )
        return False


    def add_rt_output(self, fname: str) -> None:
        """
        Add base counts from a REDItools output file.

        Parameters
        ----------
        fname : str
            Path to the REDItools output file.
        """
        self.targets.reset()
        self.exclusions.reset()
        with open_stream(fname) as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if self.do_ignore(row):
                    continue
                for nuc, count in zip(_nucs, self._counts_to_list(row[_count])):
                    key = f'{nuc}-{row[_ref]}'
                    self.counts[key] = self.counts.get(key, 0) + count

    def calc_index(self) -> dict[str, float]:
        """
        Calculate editing indices for all base transitions.

        Returns
        -------
        dict[str, float]
            A dictionary mapping transition keys (e.g., 'G-A') to editing
            indices.
        """
        indices: dict[str, float] = {}
        for idx in set(self.counts) - {f'{nuc}-{nuc}' for nuc in _nucs}:
            ref = idx[-1]
            numerator = self.counts[idx]
            denominator = self.counts.get(f'{ref}-{ref}', 0) + numerator
            if denominator == 0:
                indices[idx] = 0
            else:
                indices[idx] = 100 * numerator / denominator
        return indices

    @classmethod
    def _counts_to_list(cls, counts_str):
        pieces = counts_str[1:-1].split(', ')
        return (int(_) for _ in pieces)
