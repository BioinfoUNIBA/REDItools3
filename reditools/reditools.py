
from typing import Iterator

from reditools.alignment_manager import AlignmentManager
from reditools.compiled_position import CompiledPosition, RTResult
from reditools.compiled_reads import CompiledReads
from reditools.logger import Logger
from reditools.region import Region


class REDItools:
    """
    Main class for running REDItools analysis.

    Provides methods to set up analysis parameters and process alignment data.
    """

    def __init__(self):
        """
        Initialize REDItools with default parameters.
        """
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
        Get the current logging level.

        Returns
        -------
        str
            The current logging level.
        """
        return self._logger.level

    @log_level.setter
    def log_level(self, level: str):
        """
        Set the logging level.

        Parameters
        ----------
        level : str
            The logging level to set (e.g., 'debug', 'info', or 'silent').
        """
        self._logger = Logger(level)
        self.log = self._logger.log

    def analyze(
            self,
            alignment_manager: AlignmentManager,
            region: Region,
    ) -> Iterator[RTResult]:
        """
        Analyze a genomic region using alignment data.

        Parameters
        ----------
        alignment_manager : AlignmentManager
            The manager providing access to alignment files.
        region : Region
            The genomic region to analyze.

        Yields
        ------
        Iterator[RTResult]
            The results of the analysis for each position in the region.
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
        """
        Enable strand correction during analysis.

        Strand correction will filter reads to only those of the consensus
        strand and also report the complement of the edits and reference base
        if the strand is '-'.
        """
        self._use_strand_correction = True

    def add_reference(self, reference_fname: str):
        """
        Add a reference FASTA file for genomic reference sequences.

        Parameters
        ----------
        reference_fname : str
            The path to the reference FASTA file.
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
