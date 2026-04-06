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
        if self.stop is None:
            if self.start > 0:
                return f'{self.contig}:{self.start + 1}'
            return self.contig
        return f'{self.contig}:{self.start + 1}-{self.stop}'

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

    def contains(self, contig, position):
        """
        Determines if a given genomic location is within the region.

        Parameters:
            contig (str): Contig/Chromosome name
            position (int): Position

        Returns:
            bool
        """
        return self.contig == contig and self.start <= position < self.stop

    @staticmethod
    def from_string(region_str, alignment_file):
        contig, start, stop = Region.parse_string(region_str)
        if stop is None:
            with AlignmentFile(alignment_file, ignore_truncation=True) as bam:
                stop = bam.get_reference_length(contig)
        return Region(contig, start, stop)

    @staticmethod
    def parse_string(region_str):
        if region_str is None:
            return None
        region = re.split('[:-]', region_str)
        if not region:
            return None
        contig = region[0]
        start = 0
        stop = None

        if len(region) > 3:
            raise ValueError(f'Unrecognized format: {region_str}.')
        if len(region) > 1:
            start = Region._to_int(region[1]) - 1
            if start < 0:
                raise ValueError(
                    f'Start position ({region[1]}) must be greater than or '
                    'equal to one.',
                )
            if len(region) == 3:
                stop = Region._to_int(region[2])
                if stop <= start:
                    raise ValueError(
                        f'Stop position ({region[2]}) must be greater than or '
                        f'equal to start ({region[1]}).',
                    )
        return (contig, start, stop)

    @staticmethod
    def _to_int(number):
        if isinstance(number, str):
            return int(re.sub(r'[\s,]', '', number))
        if number is None:
            return None
        return int(number)
