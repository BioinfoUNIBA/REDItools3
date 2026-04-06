"""Wrappers for pysam files."""

from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile


class RTAlignmentFile(PysamAlignmentFile):
    """Wrapper for pysam.AlignmentFile to provide filtering on fetch."""

    def __new__(cls, *args, **kwargs):
        """
        Create a wrapper for pysam.AlignmentFile.

        Parameters:
            *args (list): Positional arguments for pysam.AlignmentFile()
            **kwargs (dict): Keyword arguments for pysam.AlignmentFile()

        Returns:
            PysamAlignmentFile
        """
        kwargs.pop('min_quality', None)
        kwargs.pop('min_length', None)
        return PysamAlignmentFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, min_quality=0, min_length=0, **kwargs):
        """
        Create a wrapper for pysam.AlignmentFile.

        This wrapper will only yield aligned reads that pass internal quality
        controls like minimum read length and minimum MAPQ.

        Parameters:
            *args (list): Positional arguments for pysam.AlignmentFile()
            min_quality (int): Minimum MAPQ
            min_length (int): Minimum read length
            **kwargs (dict): Keyword arguments for pysam.AlignmentFile()
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
        Set of read names not to be fetched.
        """
        return self._exclude_reads

    @exclude_reads.setter
    def exclude_reads(self, read_names):
        self._exclude_reads = set(read_names)
        self._checklist.append(self._check_read_name)

    def fetch(self, *args, **kwargs):
        """
        Fetch reads that pass interal quality control filters.

        Parameters:
            *args (list): Positional arguments for pysam.AlignmentFile.fetch
            *kwargs (list): Keyword arguments for pysam.AlignmentFile.fetch

        Yields:
             pysam.AlignedSegment
        """
        if 'region' in kwargs:
            kwargs['region'] = str(kwargs['region'])
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
            Lists of pysam.AlignedSegment
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

    def _check_quality(self, read):
        return read.mapping_quality >= self._min_quality

    def _check_length(self, read):
        return read.query_length >= self._min_length

    def _check_read_name(self, read):
        return read.query_name not in self._exclude_reads

    _flags_to_keep = {0, 16, 83, 99, 147, 163}

    def _check_read(self, read):
        if read.has_tag('SA'):
            return False
        if read.flag not in self._flags_to_keep:
            return False

        for check in self._checklist:
            if not check(read):
                return False
        return True
