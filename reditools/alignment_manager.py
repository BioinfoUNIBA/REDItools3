"""Wrappers for pysam files."""
from itertools import chain

from reditools.alignment_file import RTAlignmentFile


class ReadGroupIter:
    """Manages multiple fetch iterators."""

    _iter_idx = 0
    _reads_idx = 1
    _start_idx = 2

    def __init__(self, fetch_iters):
        """
        Combine multiple fetch iterators.

        Parameters:
            fetch_iters (iterable): The iterators to combine.
        """
        self._read_groups = []
        for itr in fetch_iters:
            reads = next(itr, None)
            if reads is not None:
                start = reads[0].reference_start
                self._read_groups.append({
                    self._iter_idx: itr,
                    self._reads_idx: reads,
                    self._start_idx: start,
                })

    def is_empty(self):
        """
        Check if there are still reads left.

        Returns:
            bool: True if empty, else False
        """
        return not self._read_groups

    def next(self):
        """
        Retrieve a list of reads that all start at the same position.

        Returns:
            list: Reads
        """
        position = self._find_start()
        reads = []
        for idx in range(len(self._read_groups) - 1, -1, -1):
            group = self._read_groups[idx]
            if group[self._start_idx] == position:
                reads.append(group[self._reads_idx])
                next_reads = next(group[self._iter_idx], None)
                if next_reads is None:
                    self._read_groups.pop(idx)
                else:
                    self._read_groups[idx] = {
                        self._iter_idx: group[self._iter_idx],
                        self._reads_idx: next_reads,
                        self._start_idx: next_reads[0].reference_start,
                    }
        return reads

    def _find_start(self):
        return min(group[self._start_idx] for group in self._read_groups)


class AlignmentManager:
    """
    Manage multiple RTAlignmentFiles with a single fetch.

    Attributes:
        min_quality (int): Minimum read quality (applied during add_file)
        min_length (int): Minimum read length (applied during add_file)
    """

    def __init__(self, *args, **kwargs):
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
            **self._bam_kwargs,
        )
        new_file.check_index()
        if exclude_reads:
            new_file.exclude_reads = exclude_reads
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
        rgi = ReadGroupIter(iters)
        while not rgi.is_empty():
            yield list(chain(*rgi.next()))
