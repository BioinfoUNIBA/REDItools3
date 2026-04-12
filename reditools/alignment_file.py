"""Wrappers for pysam files."""

from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile

class ReadQC:
    _flags_to_keep = {0, 16, 83, 99, 147, 163}

    def __init__(self, min_quality, min_length, exclude_set):
        self.min_quality = min_quality
        self.min_length = min_length
        self.exclude_set = exclude_set

        self.check_list = [self.check_baseline]
        if self.min_quality > 0:
            self.check_list.append(self.check_quality)
        if self.min_length > 0:
            self.check_list.append(self.check_length)
        if self.exclude_set:
            self.check_list.append(self.check_exclude_set)

    def check_baseline(self, read):
        return read.flag in self._flags_to_keep and not read.has_tag('SA')
    
    def check_quality(self, read):
        return read.mapping_quality >= self.min_quality

    def check_length(self, read):
        return read.query_length >= self.min_length

    def check_exclude_set(self, read):
        return read.query_name not in self.exclude_set

    def run_check(self, read):
        return all(_(read) for _ in self.check_list)


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
        kwargs.pop('exclude_set', None)
        return PysamAlignmentFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, min_quality=0, min_length=0, exclude_set=None, **kwargs, ):
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
        self.readqc = ReadQC(min_quality, min_length, exclude_set)

    def fetch(self, *args, **kwargs):
        """
        Fetch reads that pass interal quality control filters.

        Parameters:
            *args (list): Positional arguments for pysam.AlignmentFile.fetch
            *kwargs (list): Keyword arguments for pysam.AlignmentFile.fetch

        Yields:
             pysam.AlignedSegment
        """
        region = kwargs.get('region')
        if region is not None:
            kwargs['region'] = str(region)
        try:
            iterator = super().fetch(*args, **kwargs)
        except ValueError:
            return
        for read in iterator:
            if self.readqc.run_check(read):
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
