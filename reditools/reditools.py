"""
Analysis system for RNA editing events.

Authors:
    flat - 2017
    ahanden - 2022
"""

from reditools import utils
from reditools.compiled_reads import CompiledReads
from reditools.fasta_file import RTFastaFile
from reditools.logger import Logger
from reditools.rtchecks import RTChecks


class RTResult(object):
    """RNA editing analysis for a single base position."""

    def __init__(self, bases, strand, contig, position):
        """
        RNA editing analysis for a single base position.

        Parameters:
            bases (compiledPosition): Bases found by REDItools
            strand (str): Strand of the position
            contig (str): Contig name
            position (int): Genomic position
        """
        self.contig = contig
        self.position = position + 1
        self.bases = bases
        self.strand = strand
        self._variants = bases.get_variants()

    @property
    def variants(self):
        """
        The detected variants at this position.

        Returns:
            list
        """
        ref = self.bases.ref
        return [f'{ref}{base}' for base in self._variants]

    @property
    def mean_quality(self):
        """
        Mean read quality of the base position.

        Returns:
            int
        """
        if self.bases:
            return sum(self.bases.qualities) / len(self.bases)
        return 0

    @property
    def edit_ratio(self):
        """
        Edit ratio.

        Returns:
            float
        """
        if self._variants:
            max_edits = max(self.bases[base] for base in self._variants)
        else:
            max_edits = 0
        return max_edits / (max_edits + self.bases['REF'])

    @property
    def reference(self):
        """
        Base in the reference genome.

        Returns:
            str
        """
        return self.bases.ref

    @property
    def depth(self):
        """
        How many reads cover the position. (post filtering).

        Returns:
            int
        """
        return len(self.bases)

    @property
    def per_base_depth(self):
        """
        How many reads had each base for this position.

        Returns:
            list
        """
        return list(iter(self.bases))


class REDItools(object):
    """Analysis system for RNA editing events."""

    def __init__(self):
        """Create a new REDItools object."""
        self.hostname_string = utils.get_hostname_string()
        self._min_column_length = 1
        self._min_edits = 0
        self._min_edits_per_nucleotide = 0

        self.log_level = Logger.silent_level

        self.strand = 0
        self.usestrand_correction = False
        self.strand_confidence_threshold = 0.5

        self.min_base_quality = 30
        self.min_base_position = 0
        self.max_base_position = float('inf')

        self._rtqc = RTChecks()

        self.min_read_length = 30
        self._min_read_quality = 0

        self._target_positions = False
        self._exclude_positions = {}
        self._splice_positions = []
        self._poly_positions = []

        self.reference = None

    @property
    def poly_positions(self):
        """
        Omopolymeric positions to consider.

        Returns:
            list
        """
        return self._poly_positions

    @poly_positions.setter
    def poly_positions(self, regions):
        function = self._rtqc.check_poly_positions
        if regions:
            self._poly_positions = utils.enumerate_positions(regions)
            self._reqc.add(function)
        else:
            self._poly_positions = []
            self._rtqc.discard(function)

    @property
    def splice_positions(self):
        """
        Known splice sites.

        Returns:
            list
        """
        return self._splice_positions

    @splice_positions.setter
    def splice_positions(self, regions):
        function = self._rtqc.check_splice_positions
        if regions:
            self._splice_positions = utils.enumerate_positions(regions)
            self._rtqc.add(function)
        else:
            self._splice_positions = []
            self._rtqc.discard(function)

    @property
    def target_positions(self):
        """
        Only report results for these locations.

        Returns:
            list
        """
        return self._target_positions

    @target_positions.setter
    def target_positions(self, regions):
        if regions:
            self._target_positions = utils.enumerate_positions(regions)
        else:
            self._target_positions = False

    @property
    def exclude_positions(self):
        """
        Do not report results for these locations.

        Returns:
            list
        """
        return self.exclude_positions

    @exclude_positions.setter
    def exclude_positions(self, regions):
        self.exclude_positions = utils.enumerate_positions(regions)

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
        return self._min_read_quality  # noqa:DAR201

    @min_read_quality.setter
    def min_read_quality(self, threshold):
        self._min_read_quality = threshold
        function = self._rtqc.check_column_quality
        if self._min_read_quality > 0:
            self._rtqc.add(function)
        else:
            self._rtqc.discard(function)

    @property
    def min_column_length(self):
        """Minimum depth for a position to be reported."""
        return self._min_column_length  # noqa:DAR201

    @min_column_length.setter
    def min_column_length(self, threshold):
        self._min_column_length = threshold
        function = self._rtqc.check_column_min_length
        if threshold > 1:
            self._rtqc.add(function)
        else:
            self._rtqc.discard(function)

    @property
    def min_edits(self):
        """Minimum number of editing events for reporting."""
        return self._min_edits  # noqa:DAR201

    @min_edits.setter
    def min_edits(self, threshold):
        self._min_edits = threshold
        function = self._rtqc.check_column_edit_frequency
        if threshold > 0:
            self._rtqc.add(function)
        else:
            self._rtqc.discard(function)

    @property
    def min_edits_per_nucleotide(self):
        """Minimum number of edits for a single nucleotide for reporting."""
        return self._min_edits_per_nucleotide  # noqa:DAR201

    @min_edits_per_nucleotide.setter
    def min_edits_per_nucleotide(self, threshold):
        self._min_edits_per_nucleotide = threshold
        function = self._rtqc.check_column_min_edits
        if threshold > 0:
            self._rtqc.add(function)
        else:
            self._rtqc.discard(function)

    def analyze(self, alignment_manager, region=None):  # noqa:WPS231,WPS213
        """
        Detect RNA editing events.

        Parameters:
            alignment_manager (AlignmentManager): Source of reads
            region (Region): Where to look for edits

        Yields:
            Analysis results for each base position in region
        """
        if region is None:
            region = {}

        # Open the iterator
        self.log(
            Logger.info_level,
            'Fetching data from bams {} [REGION={}]',
            alignment_manager.file_list,
            region,
        )
        read_iter = alignment_manager.fetch_by_position(region=region)
        reads = next(read_iter, None)

        nucleotides = CompiledReads(
            self.strand,
            self.min_base_position,
            self.max_base_position,
            self.min_base_quality,
        )
        if self.reference:
            nucleotides.add_reference(self.reference)
        total = 0

        while reads is not None or not nucleotides.is_empty():
            if nucleotides.is_empty():
                self.log(
                    Logger.debug_level,
                    'Nucleotides is empty: skipping ahead',
                )
                position = alignment_manager.position
                contig = alignment_manager.contig
            else:
                position += 1

            if region.stop and position >= region.stop:
                break
            self.log(
                Logger.debug_level,
                'Analyzing position {} {}',
                contig,
                position,
            )

            # Get all the read(s) starting at position
            if reads and reads[0].reference_start == position:
                self.log(Logger.debug_level, 'Adding {} reads', len(reads))
                total += len(reads)
                nucleotides.add_reads(reads)
                reads = next(read_iter, None)

            # Process edits
            bases = nucleotides.pop(position)
            if position in self._exclude_positions.get(contig, []):
                self.log(Logger.debug_level, 'Listed exclusion - skipping')
                continue
            if self.target_positions:
                if position not in self.target_positions.get(contig, []):
                    self.log(
                        Logger.debug_level,
                        'Not listed for inclusion - skipping',
                    )
                    continue
            if bases is None:
                self.log(Logger.debug_level, 'No reads - skipping')
                continue
            column = self._get_column(position, bases, region)
            if column is None:
                self.log(Logger.debug_level, 'Bad column - skipping')
                continue
            self.log(
                Logger.debug_level,
                'Yielding output for {} reads',
                len(bases),
            )
            yield column

        self.log(
            Logger.info_level,
            '[REGION={}] {} total reads',
            region,
            total,
        )

    def use_strand_correction(self):
        """Only reports reads/positions that match `strand`."""
        self.usestrand_correction = True

    def add_reference(self, reference_fname):
        """
        Use a reference fasta file instead of reference from the BAM files.

        Parameters:
            reference_fname (str): File path to FASTA reference
        """
        self.reference = RTFastaFile(reference_fname)

    def _get_column(self, position, bases, region):
        strand = bases.get_strand(threshold=self.strand_confidence_threshold)
        if self.usestrand_correction:
            bases.filter_bystrand(strand)
        if strand == '-':
            bases.complement()

        past_start = position + 1 >= (region.start or 0)
        if past_start and bases is not None:
            if not self._rtqc.check(self, bases, region.contig, position):
                return None

        return RTResult(bases, strand, region.contig, position)


class REDItoolsDNA(REDItools):
    """
    Analysis system for editing events in DNA.

    Raises:
        ValueError: You cannot set the strand parameter using this class.
    """

    def __init__(self):
        """Create a new REDItoolsDNA object."""
        self.get_position_strand = lambda *_: '*'
        self._get_strand = lambda *_: '*'
        REDItools.__init__(self)

    def set_strand(self, strand):
        """
        Not applicable for DNA analysis.

        Parameters:
            strand (int): N/A

        Raises:
            ValueError: You cannot call this method for DNA analyses.
        """
        raise ValueError('Cannot set strand value if DNA is True')
