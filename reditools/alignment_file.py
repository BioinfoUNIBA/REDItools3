"""Wrappers for pysam files."""

from typing import Iterator

from pysam import AlignedSegment
from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile

from reditools.region import Region


class ReadQC:
    _flags_to_keep = {0, 16, 83, 99, 147, 163}

    def __init__(
            self,
            min_quality: int,
            min_length: int,
            excluded_read_names: set | None,
    ):
        self.min_quality = min_quality
        self.min_length = min_length
        self.excluded_read_names = excluded_read_names

        self.check_list = [self.check_baseline]
        if self.min_quality > 0:
            self.check_list.append(self.check_quality)
        if self.min_length > 0:
            self.check_list.append(self.check_length)
        if self.excluded_read_names:
            self.check_list.append(self.check_excluded_read_names)

    def check_baseline(self, read: AlignedSegment) -> bool:
        return read.flag in self._flags_to_keep and not read.has_tag('SA')
    
    def check_quality(self, read: AlignedSegment) -> bool:
        return read.mapping_quality >= self.min_quality

    def check_length(self, read: AlignedSegment) -> bool:
        return read.query_length >= self.min_length

    def check_excluded_read_names(self, read: AlignedSegment) -> bool:
        return read.query_name not in self.excluded_read_names  # type: ignore

    def run_check(self, read: AlignedSegment) -> bool:
        return all(_(read) for _ in self.check_list)


class RTAlignmentFile:
    """Wrapper for pysam.AlignmentFile to provide filtering on fetch."""

    def __init__(
            self,
            *args,
            min_quality: int=0,
            min_length: int=0,
            excluded_read_names: set | None=None,
            **kwargs,
    ):
        kwargs['ignore_truncation'] = True
        self.alignment_file = PysamAlignmentFile(
            *args,
            **kwargs,
        )
        self.alignment_file.check_index()
        self.readqc = ReadQC(min_quality, min_length, excluded_read_names)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.alignment_file.close()

    def fetch_by_position(
            self,
            region: Region | str,
            *args,
            **kwargs,
    ) -> Iterator[list[AlignedSegment]]:
        """
        Retrieve reads that all start at the same point on the reference.

        Parameters:
            *args (list): Positional arguments for fetch
            **kwargs (dict): Named arguments for fetch

        Yields:
            Lists of pysam.AlignedSegment
        """
        iterator = self.alignment_file.fetch(
            region=str(region, *args, **kwargs),
        )

        first_read = next(iterator, None)
        while first_read is not None and not self.readqc.run_check(first_read):
            first_read = next(iterator, None)
        if first_read is None:
            return

        reads = [first_read]
        ref_start = first_read.reference_start

        for read in iterator:
            if not self.readqc.run_check(read):
                continue
            if read.reference_start == ref_start:
                reads.append(read)
            else:
                yield reads
                reads = [read]
                ref_start = read.reference_start
        yield reads
