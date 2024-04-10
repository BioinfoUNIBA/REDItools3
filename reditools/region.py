"""Genomic Region."""

import re


class Region(object):
    """Genomic Region."""

    def __init__(self, **kwargs):
        """
        Create a new genomic region.

        Parameters:
            **kwargs (dict):
                string (str): String representation of a region
                OR
                contig (str): Contig name
                start (int): Genomic start
                stop (int): Genomic stop

        Raises:
            ValueError: The contig is missing
        """
        if 'string' in kwargs:
            region = self._parse_string(kwargs['string'])  # noqa:WPS529
            self.contig = region[0]
            self.start = region[1]
            self.stop = region[2]
        else:
            if 'contig' not in kwargs:
                raise ValueError('Region constructor requires a contig.')
            self.contig = kwargs['contig']
            self.start = self._to_int(kwargs.get('start', 1))
            self.stop = self._to_int(kwargs.get('stop', None))

    def __str__(self):
        """
        Put the region into standard string format.

        Returns:
            (str): contig:start-stop
        """
        region = self.contig
        if self.start:
            region = f'{region}:{self.start}'
            if self.stop:
                region = f'{region}-{self.stop}'
        return region

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
        if not self.stop or not self.start:
            raise IndexError('Can only split a region with a start and stop.')
        length = self.stop - self.start
        sub_regions = []
        for offset in range(0, length + 1, window):
            sub_regions.append(Region(
                contig=self.contig,
                start=self.start + offset,
                stop=self.start + offset + window,
            ))
        if self.start < length:
            sub_regions.append(Region(
                contig=self.contig,
                start=sub_regions[-1].stop,
                stop=self.stop,
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
            return None
        contig = region[0]
        start = None
        stop = None

        if len(region) > 3:
            raise ValueError(f'Unrecognized format: {region_str}.')
        if len(region) > 1:
            start = self._to_int(region[1])
            if len(region) == 3:
                stop = self._to_int(region[2])
        return (contig, start, stop)

    def _to_int(self, number):
        if isinstance(number, str):
            return int(re.sub(r'[\s,]', '', number))
        if number is None:
            return None
        return int(number)
