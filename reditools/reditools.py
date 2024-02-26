"""
Analysis system for RNA editing events.

Authors:
    flat - 2017
    ahanden - 2022
"""

import csv
from collections import defaultdict

from sortedcontainers import SortedSet

from reditools import utils
from reditools.alignment_manager import AlignmentManager
from reditools.compiled_reads import CompiledReads
from reditools.fasta_file import RTFastaFile
from reditools.file_utils import open_stream
from reditools.logger import Logger


class REDItools(object):
    """Analysis system for RNA editing events."""

    fieldnames = [
        'Region',
        'Position',
        'Reference',
        'Strand',
        'Coverage-q30',
        'MeanQ',
        'BaseCount[A,C,G,T]',
        'AllSubs',
        'Frequency',
        'gCoverage-q30',
        'gMeanQ',
        'gBaseCount[A,C,G,T]',
        'gAllSubs',
        'gFrequency',
    ]

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

        self._column_checks = set()

        self.min_read_length = 30
        self._min_read_quality = 0

        self._target_positions = False
        self._exclude_positions = {}
        self._splice_positions = []
        self._poly_positions = []

        self.reference = None

    def _check_splice_positions(self, **kwargs):
        if kwargs['position'] in self._splice_positions[kwargs['contig']]:
            self._log(
                Logger.debug_level,
                '[SPLICE_SITE] Discarding ({}, {}) because in splice site',
                kwargs['contig'],
                kwargs['position'],
            )
            return False
        return True

    def _check_poly_positions(self, **kwargs):
        contig_positions = self._omoplymeric_positions[kwargs['contig']]
        if kwargs['position'] in contig_positions:
            self._log(
                Logger.debug_level,
                '[OMOPOLYMERIC] Discarding position ({}, {})' +
                'because omopolymeric',
                kwargs['contig'],
                kwargs['position'],
            )
            return False
        return True

    def _check_column_min_length(self, bases):
        if len(bases) < self._min_column_length:
            self._log(
                Logger.debug_level,
                'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
                len(bases),
                self._min_column_length,
            )
            return False
        return True

    def _check_column_quality(self, bases):
        if bases:
            mean_q = sum(bases.qualities) / len(bases)
        else:
            mean_q = 0
        if mean_q < self._min_read_quality:
            self._log(
                Logger.debug_level,
                'DISCARD COLUMN mean_quality={}',
                mean_q,
            )
            return False
        return True

    def _check_column_edit_frequency(self, bases):
        edits_no = len(bases) - bases.base_frequency('ref')
        if edits_no < self._min_edits:
            self._log(
                Logger.debug_level,
                'DISCARDING COLUMN edits={}',
                edits_no,
            )
            return False
        return True

    def _check_columnmin_edits(self, bases):
        for num_edits in bases.getmin_edits():
            if 0 < num_edits < self._min_edits_per_nucleotide:
                self._log(
                    Logger.debug_level,
                    'DISCARDING COLUMN edits={}',
                    num_edits,
                )
                return False
        return True

    def _valid_column(self, position, bases, region):
        past_start = position + 1 >= region.start
        if past_start and bases is not None:
            return utils.check_list(self._column_checks, bases=bases)
        return False

    def load_poly_positions(self, fname):
        """
        Read omopolymeric positions from a file.

        Parameters:
            fname (str): File path
        """
        self._poly_positions = defaultdict(SortedSet)

        self._log(
            Logger.info_level,
            'Loading omopolymeric positions from file {}',
            fname,
        )

        self._log(Logger.info_level, 'Loading omopolymeric positions')

        with open_stream(fname, 'r') as stream:
            reader = csv.reader(stream, delimiter='\t')

            for fields in reader:
                if fields[0].startswith('#'):
                    continue

                contig = fields[0]
                start = int(fields[1])
                stop = int(fields[2])

                self._poly_positions[contig] |= set(range(start, stop))

        total = sum(len(pos) for pos in self._poly_positions.values())
        self._log(
            Logger.info_level,
            '{} total omopolymeric positions found.',
            total,
        )

        self._column_checks.add(self._check_poly_positions)

    def load_splicing_file(self, splicing_file, span):
        """
        Read splicing positions from a file.

        Parameters:
            splicing_file (str): File path
            span(int): Width of splice sites
        """
        self._splice_positions = defaultdict(SortedSet)
        self._log(
            Logger.info_level,
            'Loading known splice sites from file {}',
            splicing_file,
        )

        strand_map = {'-': 'D', '+': 'A'}

        with open_stream(splicing_file, 'r') as stream:
            total = 0
            total_array = defaultdict(int)
            for line in stream:
                fields = line.strip().split()

                chrom = fields[0]
                strand = fields[4]
                splice = fields[3]
                span = int(fields[1])

                total += span
                total_array[chrom] += span

                coe = -1 if strand_map.get(strand, None) == splice else 1
                new_positions = [1 + span + coe * fctr for fctr in range(span)]
                self._splice_positions[chrom] |= new_positions

        self._log(
            Logger.info_level,
            'Loaded {} positions from file {}',
            total,
            splicing_file,
        )
        self._log(Logger.info_level, 'Partial: {}', total_array)

        self._column_checks.add(self._check_splice_positions)

    def create_poly_positions(self, fname, span):
        """
        Generate omopolymeric position data.

        Parameters:
            fname (str): File path to write to
            span (int): Omopolymeric span
        """
        self._log(
            Logger.info_level,
            'Creating omopolymeric positions (span={}) from reference file',
            span,
        )

        chromosomes = self.reference.references()
        self._log(
            Logger.info_level,
            '{} chromosome names found',
            len(chromosomes),
        )

        self._log(
            Logger.info_level,
            'Writing omopolymeric positions to file: {}.',
            fname,
        )

        with open_stream(fname, 'w') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow([
                '#Chromosome',
                'Start',
                'End',
                'Length',
                'Symbol',
            ])

            for chromosome in chromosomes:
                self._log(
                    Logger.info_level,
                    'Loading reference sequence for chromosome {}',
                    chromosome,
                )
                sequence = self.reference.fetch(chromosome).lower()
                self._log(
                    Logger.info_level,
                    'Reference sequence for chromosome {} loaded (len: {})',
                    chromosome,
                    len(sequence),
                )

                equals = 0
                last = None
                for pos, base in enumerate(sequence):
                    if base == last:
                        equals += 1
                    else:
                        if equals >= span:
                            writer.writerow([
                                chromosome,
                                pos - equals,
                                pos,
                                equals,
                                last,
                            ])
                        equals = 1
                    last = base

    def load_target_positions(self, bed_file):
        """
        Read target positions from a file.

        Parameters:
            bed_file (str): Path to a BED formatted file for reading.
        """
        self._log(
            Logger.info_level,
            'Loading target positions from file {}',
            bed_file,
        )
        reader = utils.read_bed_file(bed_file)
        self._target_positions = utils.enumerate_positions(reader)

    def load_exclude_positions(self, bed_file):
        """
        Read positions to exclude from a file.

        Parameters:
            bed_file (str): Path to a BED formatted file for reading.
        """
        self._log(
            Logger.info_level,
            'Loading exclude positions from file {}',
            bed_file,
        )

        reader = utils.read_bed_file(bed_file)
        self._exclude_positions = utils.enumerate_positions(reader)

    def _next_position(self, reads, nucleotides, contig, position):
        if nucleotides.is_empty():
            self._log(
                Logger.debug_level,
                'Nucleotides is empty: skipping ahead',
            )
            position = reads[0].reference_start - 1
            contig = reads[0].reference_name
        return (contig, position + 1)

    def _get_column(self, position, bases, region):
        strand = bases.get_strand(threshold=self.strand_confidence_threshold)
        if self.usestrand_correction:
            bases.filter_bystrand(strand)
        if strand == '-':
            bases.complement()

        if not self._valid_column(position, bases, region):
            return None

        variants = bases.get_variants()
        if bases:
            mean_q = sum(bases.qualities) / len(bases)
        else:
            mean_q = 0
        return [
            position + 1,  # 1 indexed
            bases.ref,
            strand,
            len(bases),
            f'{mean_q:.2f}',
            bases.base_frequency(),
            ' '.join(sorted(variants)) if variants else '-',
            f'{bases.get_max_ratio():.2f}',
            '\t'.join(['-' for _ in range(5)]),
        ]

    def analyze(self, bam_files, region=None):
        """
        Compute editting rates for a given SAM file.

        Parameters:
            bam_files (list): Path to the SAM files
            region (dict): Must have a "contig" key and optional "start"
            and "end".

        Yields:
            Analysis results for individual base positions.
        """
        if region is None:
            region = {}

        sam_manager = AlignmentManager(
            ignore_truncation=True,
        )
        sam_manager.min_quality = self._min_read_quality
        sam_manager.min_length = self.min_read_length
        for bam in bam_files:
            sam_manager.add_file(bam)

        # Open the iterator
        self._log(
            Logger.info_level,
            'Fetching data from bams {} [REGION={}]',
            bam_files,
            region,
        )
        read_iter = sam_manager.fetch_by_position(region=region)
        reads = next(read_iter, None)

        contig = None
        position = None
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
            contig, position = self._next_position(
                reads,
                nucleotides,
                contig,
                position,
            )
            if position >= region.stop:
                break
            self._log(
                Logger.debug_level,
                'Analyzing position {} {}',
                contig,
                position,
            )

            # Get all the read(s) starting at position
            if reads is not None and reads[0].reference_start == position:
                self._log(Logger.debug_level, 'Adding {} reads', len(reads))
                total += len(reads)
                nucleotides.add_reads(reads)
                reads = next(read_iter, None)

            # Process edits
            bases = nucleotides.pop(position)
            if position in self._exclude_positions.get(contig, []):
                self._log(Logger.debug_level, 'Listed exclusion - skipping')
                continue
            if self._target_positions:
                if position not in self._target_positions.get(contig, []):
                    self._log(
                        Logger.debug_level,
                        'Not listed for inclusion - skipping',
                    )
                    continue
            if bases is None:
                self._log(Logger.debug_level, 'No reads - skipping')
                continue
            column = self._get_column(position, bases, region)
            if column is None:
                self._log(Logger.debug_level, 'Bad column - skipping')
                continue
            self._log(
                Logger.debug_level,
                'Yielding output for {} reads',
                len(bases),
            )
            yield [contig] + column

        self._log(
            Logger.info_level,
            '[REGION={}] {} total reads read',
            region,
            total,
        )

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
        self._log = self._logger.log

    @property
    def min_read_quality(self):
        """Minimum read quality for inclusion."""
        return self._min_read_quality  # noqa:DAR201

    @min_read_quality.setter
    def min_read_quality(self, threshold):
        self._min_read_quality = threshold
        if self._min_read_quality > 0:
            self._column_checks.add(self._check_column_quality)
        else:
            self._column_checks.discard(self._check_column_quality)

    @property
    def min_column_length(self):
        """Minimum depth for a position to be reported."""
        return self._min_column_length  # noqa:DAR201

    @min_column_length.setter
    def min_column_length(self, threshold):
        self._min_column_length = threshold
        if threshold > 1:
            self._column_checks.add(self._check_column_min_length)
        else:
            self._column_checks.discard(self._check_column_min_length)

    def use_strand_correction(self):
        """Only reports reads/positions that match `strand`."""
        self.usestrand_correction = True

    @property
    def min_edits(self):
        """Minimum number of editing events for reporting."""
        return self._min_edits  # noqa:DAR201

    @min_edits.setter
    def min_edits(self, threshold):
        self._min_edits = threshold
        if threshold > 0:
            self._column_checks.add(self._check_column_edit_frequency)
        else:
            self._column_checks.discard(self._check_column_edit_frequency)

    @property
    def min_edits_per_nucleotide(self):
        """Minimum number of edits for a single nucleotide for reporting."""
        return self._min_edits_per_nucleotide  # noqa:DAR201

    @min_edits_per_nucleotide.setter
    def min_edits_per_nucleotide(self, threshold):
        self._min_edits_per_nucleotide = threshold
        if threshold > 0:
            self._column_checks.add(self._check_columnmin_edits)
        else:
            self._column_checks.discard(self._check_columnmin_edits)

    def add_reference(self, reference_fname):
        """
        Use a reference fasta file instead of reference from the BAM files.

        Parameters:
            reference_fname (str): File path to FASTA reference
        """
        self.reference = RTFastaFile(reference_fname)


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
