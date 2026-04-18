"""
Analysis system for RNA editing events.

Authors:
    flat - 2017
    ahanden - 2022
"""

from reditools.compiled_position import RTResult, CompiledPosition
from reditools.compiled_reads import CompiledReads
from reditools.alignment_manager import AlignmentManager
from reditools.region import Region
from reditools.logger import Logger
from typing import Iterator


class REDItools:
    """Analysis system for RNA editing events."""

    def __init__(self):
        """Create a new REDItools object."""
        self._min_column_length = 1
        self._min_edits = 0
        self._min_edits_per_nucleotide = 0

        self._logger = Logger(Logger.silent_level)
        self.log = self._logger.log

        self.strand = 0
        self._use_strand_correction = False
        self.strand_confidence_threshold = 0.5

        self.min_base_quality = 30
        self.min_base_position = 0
        self.max_base_position = 0

        self._min_read_quality = 0

        self._specific_edits = None

        self.reference = None

    @property
    def log_level(self) -> str:
        """
        The logging level.

        Returns:
            Log level
        """
        return self._logger.level

    @log_level.setter
    def log_level(self, level: str):
        """
        Set the class logging level.

        Parameters:
            level (str): logging level
        """
        self._logger = Logger(level)
        self.log = self._logger.log

    def analyze(
            self,
            alignment_manager: AlignmentManager,
            region: Region,
    ) -> Iterator[RTResult]:
        """
        Detect RNA editing events.

        Parameters:
            alignment_manager (AlignmentManager): Source of reads
            region (Region): Where to look for edits

        Yields:
            Analysis results for each base position in region
        """
        nucleotides = CompiledReads(
            self.strand,
            self.min_base_position,
            self.max_base_position,
            self.min_base_quality,
            self.reference,
        )
        total = 0

        self.log(
            Logger.info_level,
            'Fetching reads [FILELIST={}] [REGION={}]',
            alignment_manager.file_list,
            region,
        )

        for reads in alignment_manager.fetch_by_position(region=region):
            self.log(
                Logger.debug_level,
                'Adding {} reads starting from {}:{}',
                len(reads),
                region.contig,
                reads[0].reference_start,
            )
            total += len(reads)
            nucleotides.add_reads(reads)

            next_read_start = alignment_manager.next_read_start
            if next_read_start is None or next_read_start > region.stop:
                next_read_start = region.stop

            for bases in nucleotides.pop_range(
                    reads[0].reference_start,
                    next_read_start,
            ):
                rtresult = self._process_bases(bases)
                if bases.position >= region.start:
                    self.log(
                        Logger.debug_level,
                        'Yielding output for {} reads',
                        len(rtresult),
                    )
                    yield rtresult
        self.log(
            Logger.info_level,
            '[REGION={}] {} total reads',
            region,
            total,
        )

    def use_strand_correction(self) -> None:
        """Only reports reads/positions that match `strand`."""
        self._use_strand_correction = True

    def add_reference(self, reference_fname: str):
        """
        Use a reference fasta file instead of reference from the BAM files.

        Parameters:
            reference_fname (str): File path to FASTA reference
        """
        self.reference = reference_fname

    def _process_bases(self, bases: CompiledPosition) -> RTResult:
        self.log(
            Logger.debug_level,
            'Analyzing position {} {}',
            bases.position,
            bases.contig,
        )
        if self.strand == 2:
            strand = '*'
        else:
            strand = bases.calculate_strand(
                threshold=self.strand_confidence_threshold,
            )
            if self._use_strand_correction and strand != '*':
                bases.filter_by_strand(strand)
                if strand == '-':
                    bases.complement()
        return RTResult(bases, strand)
