"""Genomic Region Collection."""
from collections import defaultdict
from typing import Iterable
from reditools.region import Region


class RegionCollection:
    """Collections of REDItools3 region objects. This class is meant to
    provide fast lookups of overlaps, and so behaves as a queue."""

    def __init__(self):
        """
        Creates a new region collection.

        """

        self._regions = defaultdict(list)
        self._index = None
        self._last_contig = None
        self._sorted = False

    def __bool__(self):
        return bool(self._regions)

    def sort(self) -> None:
        for contig, regions in self._regions.items():
            self._regions[contig] = sorted(regions)
        self._sorted = True

    def contains(self, contig: str, position: int) -> bool:
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
            self.sort()
            self._last_contig = contig
            #self._index = 0
            idx = 0
        elif contig != self._last_contig:
            self._last_contig = contig
            #self._index = 0
            idx = 0
        else:
            idx = self._index

        #for self._index, region in enumerate(
        for idx, region in enumerate(
                self._regions[contig][idx:],
                start=idx,
        ):
            if position < region.start:
                self._index = idx
                return False
            if region.start <= position < region.stop:
                self._index = idx
                return True
        self._index = len(self._regions[contig])
        return False

    def add_regions(self, regions: Iterable[Region]) -> None:
        """
        Add a list or iterable of regions to the collection.

        Parameters:
            regions (iterable): List of regions.
        """
        self._sorted = False
        for _ in regions:
            self._regions[_.contig].append(_)
