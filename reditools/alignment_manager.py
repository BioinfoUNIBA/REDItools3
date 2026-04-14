"""Wrappers for pysam files."""
from itertools import chain
from reditools.alignment_file import RTAlignmentFile

class ReadGroupIter:
    __slots__ = ('iterator', 'reads', 'reference_start')

    def __init__(self, iterator):
        self.iterator = iterator
        next(self)

    def __bool__(self):
        return bool(self.reads)

    def __next__(self):
        self.reads = next(self.iterator, None)
        if self.reads:
            self.reference_start = self.reads[0].reference_start
        else:
            self.reference_start = None
        return self.reads

class FetchGroupIter:
    """Manages multiple fetch iterators."""

    def __init__(self, fetch_iters):
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

    def __iter__(self):
        while self:
            yield next(self)

    def __bool__(self):
        return bool(self.read_groups)

    def __next__(self):
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
        return list(chain(*reads))

class AlignmentManager:
    """
    Manage multiple RTAlignmentFiles with a single fetch.

    Attributes:
        min_quality (int): Minimum read quality (applied during add_file)
        min_length (int): Minimum read length (applied during add_file)
    """

    def __init__(self, exclude_set=None, *args, **kwargs):  # noqa: WPS475
        """
        Create a new manager.

        Parameters:
            *args (list): positional arguments for pysam.AlignmentFile
                constructor
            **kwargs (dict): named arguments for pysam.AlignmentFile
                constructor
        """
        self._bam_args = args
        self._bam_kwargs = kwargs
        self._bams = []
        self.file_list = []
        self.next_read_start = None
        self.exclude_set = exclude_set

    def add_file(self, fname, exclude_reads=None):
        """
        Add an alignment file to the manager for analysis.

        Parameters:
            fname (str): Path to BAM file
            exclude_reads (set): Read names not to skip
        """
        new_file = RTAlignmentFile(
            fname,
            *self._bam_args,
            exclude_set=self.exclude_set,
            **self._bam_kwargs,
        )
        new_file.check_index()
        self._bams.append(new_file)
        self.file_list.append(fname)

    def fetch_by_position(self, *args, **kwargs):
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
