
import re
from dataclasses import dataclass

from pysam import AlignmentFile


@dataclass(slots=True, order=True, frozen=True)
class Region:
    """
    Represent a genomic region.

    Parameters
    ----------
    contig : str
        The name of the contig or chromosome.
    start : int
        The 0-based start position.
    stop : int
        The 0-based stop position (exclusive).
    """

    contig: str
    start: int
    stop: int

    def __str__(self):
        """
        Return a string representation of the region.

        Returns
        -------
        str
            The region string in 'contig:start-stop' format.
        """
        one_idx_start = self.start + 1
        if self.stop is None:
            if self.start > 0:
                return f'{self.contig}:{one_idx_start}'
            return self.contig
        return f'{self.contig}:{one_idx_start}-{self.stop}'

    def split(self, window: int) -> list['Region']:
        """
        Split the region into smaller sub-regions of a specified window size.

        Parameters
        ----------
        window : int
            The size of each sub-region.

        Returns
        -------
        list[Region]
            A list of smaller sub-regions.

        Raises
        ------
        IndexError
            If either start or stop is None.
        """
        if self.stop is None or self.start is None:
            raise IndexError('Can only split a region with a start and stop.')
        sub_regions = []
        for new_start in range(self.start, self.stop, window):
            sub_regions.append(Region(
                contig=self.contig,
                start=new_start,
                stop=min(new_start + window, self.stop)))
        return sub_regions

    @classmethod
    def from_string(cls, region_str: str, alignment_file: str) -> 'Region':
        """
        Create a Region object from a string and an alignment file.

        Parameters
        ----------
        region_str : str
            The region string in the format 'contig:start-stop',
            'contig:start', or 'contig'.
        alignment_file : str
            Path to the alignment file (BAM/CRAM) to get contig length if stop
            is missing.

        Returns
        -------
        Region
            The created Region object.

        Raises
        ------
        ValueError
            If start is less than 0 or if stop is less than or equal to start.
        """
        contig, start, stop = Region.parse_string(region_str)
        if start is None:
            start = 0
        elif start < 0:
            raise ValueError(
                f'Start position ({start}) must be greater than or '
                'equal to one.',
            )
        if stop is None:
            with AlignmentFile(alignment_file, ignore_truncation=True) as bam:
                stop = bam.get_reference_length(contig)
        if stop <= start:
            raise ValueError(
                f'Stop position ({stop}) must be greater than or '
                f'equal to start ({start}).',
            )
        return Region(contig, start, stop)

    @classmethod
    def parse_string(cls, region_str: str) -> tuple[str, int, int | None]:
        """
        Parse a region string into its components.

        Parameters
        ----------
        region_str : str
            The region string in the format 'contig:start-stop',
            'contig:start', or 'contig'.

        Returns
        -------
        tuple[str, int, int | None]
            A tuple containing (contig, start, stop). Returns None if
            region_str is None.

        Raises
        ------
        ValueError
            If the region string format is unrecognized.
        """
        if region_str is None:
            return None
        pa = re.compile(
            '(?P<contig>[^:]+)(:(?P<start>[0-9,]+)(-(?P<stop>[0-9,]+))?)?',
        )
        match = pa.fullmatch(region_str)
        if match is None:
            raise ValueError(f'Unrecognized format: {region_str}.')
        contig, start, stop = match.group('contig', 'start', 'stop')

        if start is None:
            start = 0
        else:
            start = Region._to_int(start) - 1

        if stop is not None:
            stop = Region._to_int(stop)
        return (contig, start, stop)

    @classmethod
    def _to_int(cls, number: str) -> int:
        return int(re.sub(r'[\s,]', '', number))
