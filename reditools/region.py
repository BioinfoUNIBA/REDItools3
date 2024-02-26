"""Commandline tool for REDItools."""

import re

class Region(object):
    def __init__(self, **kwargs):
        if 'string' in kwargs:
            self.contig, self.start, self.stop = self._parse_string(region_str)
        else:
            if 'contig' not in kwargs:
                raise ValueError('Region constructor requires a contig.')
            self.contig = kwargs['contig']
            if 'start' in kwargs:
                self.start = self._to_int(kwargs['start'])
            if 'stop' in kwargs:
                self.stop = self._to_int(kwargs['stop'])
        if self.start >= self.stop:
            raise IndexError(f'Start {self.start} is greater than stop {self.stop}.')

    def __str__(self):
        return f'{self.contig}:{self.start}-{self.stop}'

    def split(self, window):
        if not self.stop or not self.start:
            raise IndexError(f'Can only split a region with a start and stop.')
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
                start=self.start + offset,
                stop=self.stop,
            ))
        return sub_regions

    def _parse_string(region_str):
        """
        Parse a region string into chromosome, start, and end.

        Parameters:
            region_str (str): In the format of 'contig', 'contig:start', or 'contig:start-end'

        Returns:
            A dict of the region with keys "contig", "start", and "stop".

        Raises:
            ValueError: On improper string format
            IndexError: If the start is greater than the stop
        """
        if region_str is None:
            return None
        region = re.split('[:-]', region_str)
        if not region:
            return None
        contig = region[0]
        start = None
        stop = None

        if len(region) > 1:
            start = _to_int(region[1])
            if len(region) == 3:
                stop = _to_int(region[2])
        else:
            raise ValueError('Unrecognized format: {region_str}.')
        return (contig, start, stop)

    def _to_int(self, string):
        if type(string) == str:
            return int(re.sub(r'[\s,]', '', string))
        return int(string)
