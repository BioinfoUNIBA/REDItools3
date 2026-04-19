
from typing import Collection, Iterator

from pysam import AlignedSegment
from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile

from reditools.region import Region


class ReadQC:
    """
    Perform quality control checks on aligned reads.

    Parameters
    ----------
    min_quality : int
        Minimum mapping quality required.
    min_length : int
        Minimum read length required.
    excluded_read_names : Collection[str] | None
        A collection of read names to be excluded.
    """
    _flags_to_keep = {0, 16, 83, 99, 147, 163}

    def __init__(
            self,
            min_quality: int,
            min_length: int,
            excluded_read_names: Collection[str] | None,
    ):
        """
        Initialize the ReadQC with quality and length thresholds.

        Parameters
        ----------
        min_quality : int
            Minimum mapping quality required.
        min_length : int
            Minimum read length required.
        excluded_read_names : Collection[str] | None
            A collection of read names to be excluded.
        """
        self.min_quality = min_quality
        self.min_length = min_length
        self.excluded_read_names = excluded_read_names

        self.check_list = [self.check_baseline]
        if self.min_quality > 0:
            self.check_list.append(self.check_quality)
        if self.min_length > 0:
            self.check_list.append(self.check_length)
        if self.excluded_read_names:
            self.excluded_read_names = set(self.excluded_read_names)
            self.check_list.append(self.check_excluded_read_names)

    def check_baseline(self, read: AlignedSegment) -> bool:
        """
        Check if the read passes baseline flag and tag requirements.

        Parameters
        ----------
        read : AlignedSegment
            The aligned read to check.

        Returns
        -------
        bool
            True if the read passes, False otherwise.
        """
        return read.flag in self._flags_to_keep and not read.has_tag('SA')
    
    def check_quality(self, read: AlignedSegment) -> bool:
        """
        Check if the read passes the minimum mapping quality threshold.

        Parameters
        ----------
        read : AlignedSegment
            The aligned read to check.

        Returns
        -------
        bool
            True if the read passes, False otherwise.
        """
        return read.mapping_quality >= self.min_quality

    def check_length(self, read: AlignedSegment) -> bool:
        """
        Check if the read passes the minimum length threshold.

        Parameters
        ----------
        read : AlignedSegment
            The aligned read to check.

        Returns
        -------
        bool
            True if the read passes, False otherwise.
        """
        return read.query_length >= self.min_length

    def check_excluded_read_names(self, read: AlignedSegment) -> bool:
        """
        Check if the read name is not in the excluded list.

        Parameters
        ----------
        read : AlignedSegment
            The aligned read to check.

        Returns
        -------
        bool
            True if the read name is not excluded, False otherwise.
        """
        return read.query_name not in self.excluded_read_names  # type: ignore

    def run_check(self, read: AlignedSegment) -> bool:
        """
        Run all configured quality control checks on the read.

        Parameters
        ----------
        read : AlignedSegment
            The aligned read to check.

        Returns
        -------
        bool
            True if the read passes all checks, False otherwise.
        """
        return all(_(read) for _ in self.check_list)


class RTAlignmentFile:
    """
    A wrapper around pysam.AlignmentFile with integrated quality control.

    Parameters
    ----------
    *args
        Arguments passed to pysam.AlignmentFile.
    min_quality : int, optional
        Minimum mapping quality (default is 0).
    min_length : int, optional
        Minimum read length (default is 0).
    excluded_read_names : Collection[str] | None, optional
        A collection of read names to exclude (default is None).
    **kwargs
        Keyword arguments passed to pysam.AlignmentFile.
    """

    def __init__(
            self,
            *args,
            min_quality: int=0,
            min_length: int=0,
            excluded_read_names: Collection[str] | None=None,
            **kwargs,
    ):
        """
        Initialize the RTAlignmentFile.

        Parameters
        ----------
        *args
            Arguments passed to pysam.AlignmentFile.
        min_quality : int, optional
            Minimum mapping quality (default is 0).
        min_length : int, optional
            Minimum read length (default is 0).
        excluded_read_names : Collection[str] | None, optional
            A collection of read names to exclude (default is None).
        **kwargs
            Keyword arguments passed to pysam.AlignmentFile.
        """
        kwargs['ignore_truncation'] = True
        self.alignment_file = PysamAlignmentFile(
            *args,
            **kwargs,
        )
        self.alignment_file.check_index()
        self.readqc = ReadQC(min_quality, min_length, excluded_read_names)

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        RTAlignmentFile
            The RTAlignmentFile instance.
        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exc_type : type | None
            The exception type.
        exc_value : Exception | None
            The exception value.
        traceback : TracebackType | None
            The traceback.
        """
        self.alignment_file.close()

    def fetch_by_position(
            self,
            region: Region | str,
            *args,
            **kwargs,
    ) -> Iterator[list[AlignedSegment]]:
        """
        Fetch reads from the alignment file grouped by their reference start
        position.

        Parameters
        ----------
        region : Region | str
            The genomic region to fetch reads from.
        *args
            Arguments passed to the pysam.AlignmentFile.fetch method.
        **kwargs
            Keyword arguments passed to the pysam.AlignmentFile.fetch method.

        Yields
        ------
        list[AlignedSegment]
            A list of reads that share the same reference start position.
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
