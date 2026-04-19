from collections import defaultdict
from typing import Iterable

from reditools.region import Region


class RegionCollection:
    """
    A collection of genomic regions, providing efficient ordered lookup.
    """

    def __init__(self):
        """
        Initialize an empty RegionCollection.
        """

        self._regions = defaultdict(list)
        self._index = None
        self._last_contig = None
        self._sorted = False

    def __bool__(self):
        """
        Check whether the collection is empty.

        Returns
        -------
        bool
            True if there are regions in the collection, False otherwise.
        """
        return bool(self._regions)

    def sort(self) -> None:
        """
        Sort all regions within each contig.
        """
        for contig, regions in self._regions.items():
            self._regions[contig] = sorted(regions)
        self._sorted = True

    def contains(self, contig: str, position: int) -> bool:
        """
        Check if a given position is contained within any region of the
        collection.

        This method only works if each subsequent call is done in sorted order.
        Otherwise the output will be inconsistent.

        Parameters
        ----------
        contig : str
            The name of the contig or chromosome.
        position : int
            The 0-based position to check.

        Returns
        -------
        bool
            True if the position is within a region, False otherwise.
        """
        if not self._sorted:
            self.sort()
            self._last_contig = contig
            idx = 0
        elif contig != self._last_contig:
            self._last_contig = contig
            idx = 0
        else:
            idx = self._index


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
        self._index = len(self.get_contig(contig))
        return False

    def add_regions(self, regions: Iterable[Region]) -> None:
        """
        Add multiple regions to the collection.

        Parameters
        ----------
        regions : Iterable[Region]
            An iterable of Region objects to add.
        """
        self._sorted = False
        for _ in regions:
            self._regions[_.contig].append(_)

    def get_contig(self, contig: str) -> list[Region]:
        """
        Retrieve the regions for a specific chromosome/contig.

        Parameters
        ----------
        contig : str
            Chromosome/contig to retrieve.

        Returns
        -------
        list[Region]
            Regions found in the given contig.
        """
        return self._regions[contig]
