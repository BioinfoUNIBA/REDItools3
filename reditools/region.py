"""Genomic Region."""

import re


class Region(object):
    """Genomic Region."""

    def __init__(self, **kwargs):
        """
        Create a new genomic region.

        Parameters:
            **kwargs (dict):
                string (str): String representation of a region in
                              samtools-compatible notation
                OR
                contig (str): Contig name
                start (int): Genomic start (zero index, inclusive)
                stop (int): Genomic stop (zero index, exclusive)

        Raises:
            ValueError: The contig is missing
        """
        if 'string' in kwargs:
            region = self._parse_string(kwargs['string'])  # noqa:WPS529
            self.contig, self.start, self.stop = region
        else:
            if 'contig' not in kwargs:
                raise ValueError('Region constructor requires a contig.')
            self.contig = kwargs['contig']
            self.start = self._to_int(kwargs.get('start', 0))
            self.stop = self._to_int(kwargs.get('stop', None))

    def __str__(self):
        """
        Put the region into standard string format.

        Returns:
            (str): contig:start-stop
        """
        if self.start >= 0:
            if self.stop is not None:
                return f'{self.contig}:{self.start + 1}-{self.stop + 1}'
            return f'{self.contig}:{self.start + 1}'
        return self.contig

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
        length = self.stop - self.start
        sub_regions = []
        for offset in range(0, length + 1, window):
            sub_regions.append(Region(
                contig=self.contig,
                start=self.start + offset,
                stop=min(self.start + offset + window, self.stop),
            ))
        return sub_regions

    def enumerate(self):
        """
        Convert a list of regions into a list of individual positions.

        Returns:
            Set enumerating the individual positions.
        """
        return set(range(self.start, self.stop))

    def contains(self, contig, position):
        """
        Determines if a given genomic location is within the region.

        Parameters:
            contig (str): Contig/Chromosome name
            position (int): Position

        Returns:
            bool
        """
        if self.contig != contig:
            return False
        left = self.start is None or self.start <= position
        right = self.stop is None or position < self.stop
        return left and right

    def _parse_string(self, region_str):
        if region_str is None:
            return None
        region = re.split('[:-]', region_str)
        if not region:
            return [None, None, None]
        contig = region[0]
        start = 0
        stop = None

        if len(region) > 3:
            raise ValueError(f'Unrecognized format: {region_str}.')
        if len(region) > 1:
            start = self._to_int(region[1]) - 1
            if len(region) == 3:
                stop = self._to_int(region[2]) - 1
        return (contig, start, stop)

    def _to_int(self, number):
        if isinstance(number, str):
            return int(re.sub(r'[\s,]', '', number))
        if number is None:
            return None
        return int(number)
