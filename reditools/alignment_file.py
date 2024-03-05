"""Wrappers for pysam files."""

from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile


class RTAlignmentFile(PysamAlignmentFile):
    """Wrapper for pysam.AlignmentFile to provide filtering on fetch."""

    def __new__(cls, *args, **kwargs):
        """
        Create a wrapper for pysam.AlignmentFile.

        Parameters:
            *args (list): Positional arguments for pysam.FastaFile()
            **kwargs (dict): Keyword arguments for pysam.FastaFile()

        Returns:
            PysamAlignmentFile
        """
        kwargs.pop('min_quality', None)
        kwargs.pop('min_length', None)
        return PysamAlignmentFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, min_quality=0, min_length=0, **kwargs):
        """
        Create a wrapper for pysam.AlignmentFile.

        Parameters:
            *args (list): Positional arguments for pysam.FastaFile()
            min_quality (int): Minimum read quality
            min_length (int): Minimum read length
            **kwargs (dict): Keyword arguments for pysam.FastaFile()
        """
        PysamAlignmentFile.__init__(self)

        self._checklist = []

        if min_quality > 0:
            self._min_quality = min_quality
            self._checklist.append(self._check_quality)

        if min_length > 0:
            self._min_length = min_length
            self._checklist.append(self._check_length)

    @property
    def exclude_reads(self):
        """
        Names of reads not to be fetched.

        Returns:
            iterable
        """
        return self._exclude_reads

    @exclude_reads.setter
    def exclude_reads(self, read_names):
        """
        Provide a list of read names to be skipped during fetch.

        Parameters:
            read_names (iterable): Reads to skip
        """
        self._exclude_reads = set(read_names)
        self._checklist.append(self._check_read_name)

    def fetch(self, *args, **kwargs):
        """
        Fetch reads aligned in a region.

        Parameters:
            *args (list): Positional arguments for pysam.FastaFile.fetch
            *kwargs (list): Keyword arguments for pysam.FastaFile.fetch

        Yields:
            Reads
        """
        if 'region' in kwargs:
            kwargs['region'] = str(kwargs['region'])  # noqa:WPS529
        try:
            iterator = super().fetch(*args, **kwargs)
        except ValueError:
            return
        for read in iterator:
            if self._check_read(read):
                yield read

    def fetch_by_position(self, *args, **kwargs):
        """
        Retrieve reads that all start at the same point on the reference.

        Parameters:
            *args (list): Positional arguments for fetch
            **kwargs (dict): Named arguments for fetch

        Yields:
            Lists containing reads
        """
        iterator = self.fetch(*args, **kwargs)

        first_read = next(iterator, None)
        if first_read is None:
            return

        reads = [first_read]
        ref_start = first_read.reference_start

        for read in iterator:
            if read.reference_start == ref_start:
                reads.append(read)
            else:
                yield reads
                reads = [read]
                ref_start = read.reference_start
        yield reads

    # 77: NOT_MAPPED
    # 141: NOT_MAPPED
    # 512: QC_FAIL
    # 256: IS_SECONDARY
    # 2048: IS_SUPPLEMENTARY
    # 1024: IS_DUPLICATE
    _flags_to_toss = {77, 141, 512, 256, 2048, 1024}
    _paired_flags_to_keep = {99, 147, 83, 163}

    def _check_quality(self, read):
        return read.mapping_quality >= self._min_quality

    def _check_length(self, read):
        return read.query_length >= self._min_length

    def _check_read_name(self, read):
        return read.query_name not in self._exclude_reads

    def _check_read(self, read):
        if read.has_tag('SA'):
            return False
        if read.flag in self._flags_to_toss:
            return False
        if read.is_paired and read.flag not in self._paired_flags_to_keep:
            return False

        for check in self._checklist:
            if not check(read):
                return False
        return True
