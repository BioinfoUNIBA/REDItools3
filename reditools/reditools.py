"""
Analysis system for RNA editing events.

Authors:
    flat - 2017
    ahanden - 2022
"""

from reditools.compiled_reads import CompiledReads
from reditools.fasta_file import RTFastaFile
from reditools.logger import Logger
from reditools.region_collection import RegionCollection
from reditools import rtchecks


class REDItools:
    """Analysis system for RNA editing events."""

    def __init__(self):
        """Create a new REDItools object."""
        self._min_column_length = 1
        self._min_edits = 0
        self._min_edits_per_nucleotide = 0

        self.log_level = Logger.silent_level

        self.strand = 0
        self._use_strand_correction = False
        self.strand_confidence_threshold = 0.5

        self.min_base_quality = 30
        self.min_base_position = 0
        self.max_base_position = 0

        self._rtqc = rtchecks.RTChecks()

        self._min_read_quality = 0

        self._target_regions = RegionCollection()
        self._exclude_regions = RegionCollection()
        self._specific_edits = None

        self.reference = None

        self._include_refs = None

    @property
    def include_refs(self):
        """
        Genome reference bases to report on.

        Returns:
            list
        """
        return self._include_refs

    @property
    def specific_edits(self):
        """
        Specific edit events to report.

        Returns:
            set
        """
        return self._specific_edits

    @specific_edits.setter
    def specific_edits(self, edits):
        if list(edits) == ["ALL"]:
            self._specific_edits = set()
        else:
            self._specific_edits = set(edits)

        for alt in self._specific_edits:
            if not self._verify_alt(alt):
                raise Exception(
                        f'Specific edit "{alt}" is not valid. ' +
                        'Edits must be two character strings of ATCG.')

        qc_check = rtchecks.check_specific_alts
        if self._specific_edits:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    def _verify_alt(self, alt):
        if not isinstance(alt, str):
            return False
        if len(alt) != 2:
            return False
        if alt[0] not in 'ATCG' and alt[1] not in 'ATCG':
            return False
        return True

    @property
    def target_regions(self):
        """
        Only report results for these locations.

        Returns:
            list
        """
        return self._target_regions

    def add_target_regions(self, regions):
        """
        Only report results for these locations.

        Parameters:
            regions (iterable): List of Region objects.
        """
        if regions:
            self._target_regions.add_regions(regions)
            self._rtqc.add(rtchecks.check_target_positions)

    @property
    def log_level(self):
        """
        The logging level.

        Returns:
            Log level
        """
        return self._log_level

    @log_level.setter
    def log_level(self, level):
        """
        Set the class logging level.

        Parameters:
            level (str): logging level
        """
        self._logger = Logger(level)
        self.log = self._logger.log

    @property
    def min_read_quality(self):
        """Minimum read quality for inclusion."""
        return self._min_read_quality

    @min_read_quality.setter
    def min_read_quality(self, threshold):
        self._min_read_quality = threshold
        qc_check = rtchecks.check_column_quality
        if self._min_read_quality > 0:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    @property
    def min_column_length(self):
        """Minimum depth for a position to be reported."""
        return self._min_column_length

    @min_column_length.setter
    def min_column_length(self, threshold):
        self._min_column_length = threshold
        qc_check = rtchecks.check_column_min_length
        if threshold > 1:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    @property
    def min_edits(self):
        """Minimum number of editing events for reporting."""
        return self._min_edits

    @min_edits.setter
    def min_edits(self, threshold):
        self._min_edits = threshold
        qc_check = rtchecks.check_column_edit_frequency
        if threshold > 0:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    @property
    def min_edits_per_nucleotide(self):
        """Minimum number of edits for a single nucleotide for reporting."""
        return self._min_edits_per_nucleotide

    @min_edits_per_nucleotide.setter
    def min_edits_per_nucleotide(self, threshold):
        self._min_edits_per_nucleotide = threshold
        qc_check = rtchecks.check_column_min_edits
        if threshold > 0:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    @property
    def exclude_regions(self):
        """Regions to exclude from analysis"""
        return self._exclude_regions

    def add_exclude_regions(self, regions):
        """
        Regions to exclude from analysis

        Parameters:
            regions (iterable): List of Region objects.
        """
        if regions:
            self._exclude_regions.add_regions(regions)
            self._rtqc.add(rtchecks.check_exclusions)

    @property
    def max_alts(self):
        """Maximum number of alternative bases for a position."""
        return self._max_alts

    @max_alts.setter
    def max_alts(self, max_alts):
        self._max_alts = max_alts
        qc_check = rtchecks.check_max_alts
        if max_alts < 3:
            self._rtqc.add(qc_check)
        else:
            self._rtqc.discard(qc_check)

    def analyze(self, alignment_manager, region):
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
        )
        if self.reference:
            nucleotides.add_reference(self.reference)
        total = 0

        self.log(
            Logger.info_level,
            'Fetching reads [FILELIST={}] [REGION={}]',
            alignment_manager.file_list,
            region,
        )
        read_iter = alignment_manager.fetch_by_position(region=region)

        reads = next(read_iter, None)
        while reads is not None:
            position = reads[0].reference_start
            self.log(
                Logger.debug_level,
                'Adding {} reads starting from {}:{}',
                len(reads),
                region.contig,
                position,
            )
            total += len(reads)
            nucleotides.add_reads(reads)

            reads = next(read_iter, None)
            if reads is None or reads[0].reference_start > region.stop:
                next_read_start = region.stop
            else:
                next_read_start = reads[0].reference_start

            for position in range(position, next_read_start):
                if nucleotides.is_empty():
                    break
                bases = nucleotides.pop(position)
                if position < region.start:
                    continue
                self.log(
                    Logger.debug_level,
                    'Analyzing position {} {}',
                    region.contig,
                    position,
                )
                if bases is None:
                    self.log(Logger.debug_level, 'DISCARD COLUMN no reads')
                    continue
                bases.calculate_strand(
                    threshold=self.strand_confidence_threshold,
                )
                if self._use_strand_correction:
                    bases.filter_by_strand()
                    if bases.strand == '-':
                        bases.complement()
                if self._rtqc.check(self, bases):
                    self.log(
                        Logger.debug_level,
                        'Yielding output for {} reads',
                        len(bases),
                    )
                    yield bases
        self.log(
            Logger.info_level,
            '[REGION={}] {} total reads',
            region,
            total,
        )

    def use_strand_correction(self):
        """Only reports reads/positions that match `strand`."""
        self._use_strand_correction = True

    def add_reference(self, reference_fname):
        """
        Use a reference fasta file instead of reference from the BAM files.

        Parameters:
            reference_fname (str): File path to FASTA reference
        """
        self.reference = RTFastaFile(reference_fname)
