"""Wrappers for pysam files."""
from itertools import chain
from typing import Collection, Iterable, Iterator

from pysam import AlignedSegment

from reditools.alignment_file import RTAlignmentFile


class ReadGroupIter:
    __slots__ = ('iterator', 'reads', 'reference_start')

    def __init__(self, iterator: Iterator):
        self.iterator = iterator
        next(self)

    def __bool__(self):
        return bool(self.reads)

    def __next__(self) -> list[AlignedSegment] | None:
        self.reads = next(self.iterator, None)
        if self.reads:
            self.reference_start = self.reads[0].reference_start
        else:
            self.reference_start = None
        return self.reads

class FetchGroupIter:
    """Manages multiple fetch iterators."""

    def __init__(self, fetch_iters: list[Iterator]):
        """
        Combine multiple fetch iterators.

        Parameters:
            fetch_iters (iterable): The iterators to combine.
        """
        self.read_groups = []
        for iterator in fetch_iters:
            rgi = ReadGroupIter(iterator)
            if rgi:
                self.read_groups.append(rgi)

    def __iter__(self) -> Iterator[list[AlignedSegment]]:
        while self:
            yield next(self)

    def __bool__(self):
        return bool(self.read_groups)

    def __next__(self) -> list[AlignedSegment]:
        """
        Retrieve a list of reads that all start at the same position.

        Returns:
            list: Reads
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
    def __init__(
            self,
            excluded_read_names: Collection[str] | None=None,
            min_quality: int=0,
            min_length: int=0,
    ):  # noqa: WPS475
        self._bams: list[RTAlignmentFile] = []
        self.file_list: list[str] = []
        self.next_read_start: int | None = None
        self.excluded_read_names = excluded_read_names
        self.min_quality = min_quality
        self.min_length = min_length

    def add_file(self, fname: str):
        """
        Add an alignment file to the manager for analysis.

        Parameters:
            fname (str): Path to BAM file
            exclude_reads (set): Read names not to skip
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
        Perform combine fetch_by_position for all managed files.

        Parameters:
            *args (list): Positional arguments for
                RTAlignmentFile.fetch_by_position
            **kwargs (dict): Named arguments for
                RTAlignmentFile.fetch_by_position

        Yields:
            list: reads from all managed files that begin at the same position.
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
