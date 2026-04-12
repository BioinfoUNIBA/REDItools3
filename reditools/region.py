"""Genomic Region."""

import re
from dataclasses import dataclass

from pysam import AlignmentFile


@dataclass(slots=True, order=True)
class Region:
    contig: str
    start: int
    stop: int
    """Genomic Region."""

    def __str__(self):
        """
        Put the region into standard string format.

        Returns:
            (str): contig:start-stop
        """
        one_idx_start = self.start + 1
        if self.stop is None:
            if self.start > 0:
                return f'{self.contig}:{one_idx_start}'
            return self.contig
        return f'{self.contig}:{one_idx_start}-{self.stop}'

    def split(self, window):
        """
        Split the region into a list of smaller regions.

        Parameters:
            window (int): The size of the sub regions in bp

        Returns:
            list

        Raises:
            IndexError: The region is missing a start or stop
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
    def from_string(cls, region_str, alignment_file):
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
    def parse_string(cls, region_str):
        if region_str is None:
            return None
        pa = re.compile(
            '(?P<contig>[^:]+)(:(?P<start>[0-9,]+)(-(?P<stop>[0-9,]+))?)?',
        )
        match = pa.fullmatch(region_str)

        try:
            contig, start, stop = match.group('contig', 'start', 'stop')
        except AttributeError as exc:
            raise ValueError(f'Unrecognized format: {region_str}.') from exc

        if start is None:
            start = 0
        else:
            start = Region._to_int(start) - 1

        if stop is not None:
            stop = Region._to_int(stop)
        return (contig, start, stop)

    @classmethod
    def _to_int(cls, number):
        if isinstance(number, str):
            return int(re.sub(r'[\s,]', '', number))
        if number is None:
            return None
        return int(number)
