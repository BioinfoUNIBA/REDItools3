from itertools import chain
from typing import Collection, Iterable, Iterator

from pysam import AlignedSegment

from reditools.alignment_file import RTAlignmentFile


class ReadGroupIter:
    """
    Iterator over groups of reads sharing the same reference start position.

    Parameters
    ----------
    iterator : Iterator
        An iterator yielding lists of AlignedSegment objects.
    """
    __slots__ = ('iterator', 'reads', 'reference_start')

    def __init__(self, iterator: Iterator):
        """
        Initialize the ReadGroupIter.

        Parameters
        ----------
        iterator : Iterator
            An iterator yielding lists of AlignedSegment objects.
        """
        self.iterator = iterator
        next(self)

    def __bool__(self):
        """
        Check if there are more reads.

        Returns
        -------
        bool
            True if there are more reads, False otherwise.
        """
        return bool(self.reads)

    def __next__(self) -> list[AlignedSegment] | None:
        """
        Get the next group of reads.

        Returns
        -------
        list[AlignedSegment] | None
            The next list of reads, or None if the iterator is exhausted.
        """
        self.reads = next(self.iterator, None)
        if self.reads:
            self.reference_start = self.reads[0].reference_start
        else:
            self.reference_start = None
        return self.reads

class FetchGroupIter:
    """
    Iterator that merges multiple ReadGroupIter objects, yielding reads grouped
    by position.

    Parameters
    ----------
    fetch_iters : list[Iterator]
        A list of iterators, each yielding reads from an alignment file.
    """

    def __init__(self, fetch_iters: list[Iterator]):
        """
        Initialize the FetchGroupIter.

        Parameters
        ----------
        fetch_iters : list[Iterator]
            A list of iterators, each yielding reads from an alignment file.
        """
        self.read_groups = []
        for iterator in fetch_iters:
            rgi = ReadGroupIter(iterator)
            if rgi:
                self.read_groups.append(rgi)

    def __iter__(self) -> Iterator[list[AlignedSegment]]:
        """
        Return the iterator object itself.

        Returns
        -------
        Iterator[list[AlignedSegment]]
            The iterator object.
        """
        while self:
            yield next(self)

    def __bool__(self):
        """
        Check if there are more read groups.

        Returns
        -------
        bool
            True if there are more read groups, False otherwise.
        """
        return bool(self.read_groups)

    def __next__(self) -> list[AlignedSegment]:
        """
        Get the next group of reads from all alignment files for the same
        position.

        Returns
        -------
        list[AlignedSegment]
            A concatenated list of reads from all files at the current minimum
            position.
        """
        position = min(_.reference_start for _ in self.read_groups)
        reads = []
        for idx, rgi in reversed(list(enumerate(self.read_groups))):
            if rgi.reference_start != position:
                continue
            reads.append(rgi.reads)
            if next(rgi) is None:
                self.read_groups.pop(idx)
        return list(chain(*reads))  # type: ignore

class AlignmentManager:
    """
    Manage multiple alignment files and provide unified access to reads by
    position.

    Parameters
    ----------
    excluded_read_names : Collection[str] | None, optional
        Read names to exclude from analysis (default is None).
    min_quality : int, optional
        Minimum mapping quality (default is 0).
    min_length : int, optional
        Minimum read length (default is 0).
    """
    def __init__(
            self,
            excluded_read_names: Collection[str] | None=None,
            min_quality: int=0,
            min_length: int=0,
    ):  # noqa: WPS475
        """
        Initialize the AlignmentManager.

        Parameters
        ----------
        excluded_read_names : Collection[str] | None, optional
            Read names to exclude from analysis (default is None).
        min_quality : int, optional
            Minimum mapping quality (default is 0).
        min_length : int, optional
            Minimum read length (default is 0).
        """
        self._bams: list[RTAlignmentFile] = []
        self.file_list: list[str] = []
        self.next_read_start: int | None = None
        self.excluded_read_names = excluded_read_names
        self.min_quality = min_quality
        self.min_length = min_length

    def add_file(self, fname: str):
        """
        Add an alignment file to the manager.

        Parameters
        ----------
        fname : str
            Path to the alignment file (BAM).
        """
        new_file = RTAlignmentFile(
            fname,
            excluded_read_names=self.excluded_read_names,
            min_length=self.min_length,
            min_quality=self.min_quality,
        )
        self._bams.append(new_file)
        self.file_list.append(fname)

    def fetch_by_position(
            self,
            *args,
            **kwargs,
    ) -> Iterable[list[AlignedSegment]]:
        """
        Fetch reads from all managed files, grouped by position.

        Parameters
        ----------
        *args
            Arguments passed to the fetch method of pysam.AlignmentFile.
        **kwargs
            Keyword arguments passed to the fetch method of
            pysam.AlignmentFile.

        Yields
        ------
        list[AlignedSegment]
            A list of reads from all files at each genomic position.
        """
        iters = [bam.fetch_by_position(*args, **kwargs) for bam in self._bams]
        fgi = FetchGroupIter(iters)
        if not fgi:
            return
        read_group = next(fgi)
        self.next_read_start = read_group[0].reference_start
        for next_read_group in fgi:
            self.next_read_start = next_read_group[0].reference_start
            yield read_group
            read_group = next_read_group
        self.next_read_start = None
        yield read_group
