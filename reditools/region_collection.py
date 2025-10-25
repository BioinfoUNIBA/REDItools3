"""Genomic Region."""

import re

from reditools.region import Region

class RegionCollection(object):
    """Collections of REDItools3 region objects. This class is meant to
    provide fast lookups of overlaps, and so behaves as a queue."""

    def __init__(self):
        """
        Creates a new region collection.

        """

        self._regions = {}
        self._last_index = None
        self._last_contig = None
        self._sorted = False

    def _sort(self):
        for contig, regions in self._regions.items():
            self._regions[contig] = sorted(regions, key=lambda _: (_.start, _.stop))
        self._sorted = True

    def contains(self, contig, position):
        """
        Checks whether the given position or range overlaps with the
        collection.

        Parameters:
            contig (str): Chromomsome/contig name from reference.
            position_start (int): Position to check for.

        Returns:
            True if there is an overlap, False otherwise.
        """
        if not self._sorted:
            self._sort()

        if contig != self._last_contig:
            self._last_contig = contig
            self._last_index = 0

        for i in range(self._last_index, len(self._regions.get(contig, []))):
            self._last_index = i
            region = self._regions[contig][i]
            if position < region.start:
                return False
            if position >= region.start and position < region.stop:
                return True
        self._last_index += 1

        return False

    def addRegion(self, region):
        """
        Add a region to the collection.
        
        Parameters:
            region (Region): region to add.
        """
        self._sorted = False
        if region.contig not in self._regions:
            self._regions[region.contig] = []
        self._regions[region.contig].append(region)

    def addRegions(self, regions):
        """
        Add a list or iterable of regions to the collection.

        Parameters:
            regions (iterable): List of regions.
        """
        for r in regions:
            self.addRegion(r)
