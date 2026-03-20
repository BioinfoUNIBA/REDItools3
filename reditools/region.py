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
            if 'stop' in kwargs:
                self.stop = self._to_int(kwargs['stop'])
            else:
                self.stop = None

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
        return f'{self.contig}:{self.start + 1}-{self.stop + 1}'

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
        start = 0
        stop = None

        if len(region) > 3:
            raise ValueError(f'Unrecognized format: {region_str}.')
        if len(region) > 1:
            start = self._to_int(region[1]) - 1
            if start < 0:
                raise ValueError(
                    f'Start position ({region[1]}) must be greater than or '
                    'equal to one.',
                )
            if len(region) == 3:
                stop = self._to_int(region[2]) - 1
                if stop <= start:
                    raise ValueError(
                        f'Stop position ({region[2]}) must be greater than or '
                        f'equal to start ({region[1]}).',
                    )
        return (contig, start, stop)

    def _to_int(self, number):
        if isinstance(number, str):
            return int(re.sub(r'[\s,]', '', number))
        if number is None:
            return None
        return int(number)
