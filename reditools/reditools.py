#!/usr/bin/env python
'''
Analysis system for RNA editing events.
Authors:
    flat - 2017
    ahanden - 2022
'''

from __future__ import division

from collections import defaultdict
from datetime import datetime
import csv
import sys

from sortedcontainers import SortedSet

from reditools import utils
from reditools.compiled_position import CompiledReads
from reditools.alignment import FastaFile, AlignmentManager

# pylint: disable=too-many-instance-attributes
class REDItools:
    '''
    Analysis system for RNA editing events.
    '''

    _DEFAULT_BASE_QUALITY = 30

    _complement_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

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
        'gFrequency']

    @staticmethod
    def _default(response):
        '''
        Returns a function reference to a method that always returns the same response.

        Parameters:
            response (mixed): Whatever the function should return

        Returns:
            A function reference.
        '''
        return lambda *x: response

    def _log_all(self, level, message, *args):
        '''
        Writes messages to stderr.

        Parameters:
            level (str): The log level
            message (str): The message text
            *args: Additional argument to pass to message.format
        '''
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        message = message.format(*args)
        sys.stderr.write(f'{timestamp} [{self._hostname_string}] [{level}] {message}\n')

    def _log_verbose(self, level, message, *args):
        '''
        Writes INFO level messages to stderr.

        Parameters:
            level (str): The log level. If level is not "INFO", nothing will be written
            message (str): The message text
            *args: Additional argument to pass to message.format
        '''

        if level == "INFO":
            self._log_all(level, message, *args)

    def _check_splice_positions(self, **kwargs):
        '''
        Determines if the given base is in a known splice site.

        Parameters:
            position (int): Base site to check
            contig (str): Chromsome

        Returns:
            True if the location is in a splice site, else False
        '''
        if kwargs["position"] in self._splice_positions[kwargs["contig"]]:
            self._log('DEBUG', '[SPLICE_SITE] Discarding ({}, {}) because in splice site',
                      kwargs["contig"], kwargs["position"])
            return False
        return True

    def _check_omopolymeric_positions(self, **kwargs):
        '''
        Determines if the given base is in a omopolymeric site.

        Parameters:
            position (int): Base site to check
            contig (str): Chromsome

        Returns:
            True if the location is in a omopolymeric site, else False
        '''

        if kwargs["position"] in self._omopolymeric_positions[kwargs["contig"]]:
            self._log('DEBUG', '[OMOPOLYMERIC] Discarding position ({}, {}) because omopolymeric',
                      kwargs["contig"], kwargs["position"])
            return False
        return True

    def _check_column_min_length(self, bases):
        if len(bases) < self._min_column_length:
            self._log('DEBUG',
                      'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
                      len(bases),
                      self._min_column_length)
            return False
        return True

    def _check_column_quality(self, bases):
        mean_q = bases.mean_quality()
        if mean_q < self._min_read_quality:
            self._log('DEBUG', 'DISCARD COLUMN mean_quality={}', mean_q)
            return False
        return True

    def _check_column_edit_frequency(self, bases):
        edits_no = len(bases) - bases.ref_frequency()
        if edits_no < self._min_edits:
            self._log('DEBUG', 'DISCARDING COLUMN edits={}', edits_no)
            return False
        return True

    def _check_column_min_edits(self, bases):
        for i in bases.get_min_edits():
            if 0 < i < self._min_edits_per_nucleotide:
                self._log('DEBUG',
                          'DISCARDING COLUMN edits={}',
                          i)
                return False
        return True

    def _valid_column(self, position, bases, region):
        return position + 1 >= region.get("start", 0) and \
               bases is not None and \
               utils.check_list(self._column_checks, bases=bases)

    def load_omopolymeric_positions(self, fname):
        '''
        Reads omopolymeric positions from a file.
        '''
        self._omopolymeric_positions = defaultdict(SortedSet)

        self._log('INFO',
                  'Loading omopolymeric positions from file {}',
                  fname)

        self._log('INFO', 'Loading omopolymeric positions')

        with utils.open_stream(fname, 'r') as stream:
            reader = csv.reader(stream, delimiter="\t")

            for fields in reader:
                if fields[0].startswith('#'):
                    continue

                contig = fields[0]
                start = int(fields[1])
                stop = int(fields[2])

                self._omopolymeric_positions[contig] |= [*range(start, stop)]

        total = sum(len(i) for i in self._omopolymeric_positions.values())
        self._log('INFO', '{} total omopolymeric positions found.', total)

        self._column_checks.add(self._check_omopolymeric_positions)

    def load_splicing_file(self, splicing_file, splicing_span):
        '''
        Reads splicing positions from a file.
        '''

        self._splice_positions = defaultdict(SortedSet)
        self._log('INFO',
                  'Loading known splice sites from file {}',
                  splicing_file)

        strand_map = {'-': 'D', '+': 'A'}

        with utils.open_stream(splicing_file, "r") as stream:
            total = 0
            total_array = defaultdict(int)
            for line in stream:
                fields = line.strip().split()

                chrom, strand, splice, span = fields[0], fields[4], fields[3], int(fields[1])

                total += splicing_span
                total_array[chrom] += splicing_span

                coe = -1 if strand_map.get(strand, None) == splice else 1
                self._splice_positions[chrom] |= [1 + span + coe * j \
                    for j in range(splicing_span)]

        self._log('INFO', 'Loaded {} positions from file {}', total, splicing_file)
        self._log('INFO', 'Partial: {}', total_array)

        self._column_checks.add(self._check_splice_positions)

    def create_omopolymeric_positions(self, fname, span):
        '''
        Generates omopolymeric position data.

        Parameters:
            omopolymeric_file (str): File path to write to
            span (int): Omopolymeric span
        '''
        self._log('INFO',
                  'Creating omopolymeric positions (span={}) from reference file',
                  span)

        chromosomes = self.reference.references()
        self._log('INFO', '{} chromosome names found', len(chromosomes))

        self._log('INFO', 'Writing omopolymeric positions to file: {}.', fname)

        with utils.open_stream(fname, 'w') as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "#Chromosome",
                "Start",
                "End",
                "Length",
                "Symbol"])

            for chromosome in chromosomes:
                self._log('INFO',
                          'Loading reference sequence for chromosome {}',
                          chromosome)
                sequence = self.reference.fetch(chromosome).lower()
                self._log('INFO',
                          'Reference sequence for chromosome {} loaded (len: {})',
                          chromosome, len(sequence))

                equals = 0
                last = None
                for pos, base in enumerate(sequence):
                    if base == last:
                        equals += 1
                    else:
                        if equals >= span:
                            writer.writerow([chromosome, pos - equals, pos, equals, last])
                        equals = 1
                    last = base

    def load_target_positions(self, bed_file):
        '''
        Reads target positions from a file.

        Parameters:
            bed_file (str): Path to a BED formatted file for reading.
        '''
        self._log('INFO',
                  'Loading target positions from file {}',
                  bed_file)
        reader = utils.read_bed_file(bed_file)
        self._target_positions = utils.enumerate_positions(reader)

    def load_exclude_positions(self, bed_file):
        '''
        Reads positions to exclude from a file.

        Parameters:
            bed_file (str) Path to a BED formatted file for reading.
        '''
        self._log('INFO',
                  'Loading exclude positions from file {}',
                  bed_file)

        reader = utils.read_bed_file(bed_file)
        self._exclude_positions = utils.enumerate_positions(reader)

    def _next_position(self, reads, nucleotides, contig, position):
        if nucleotides.is_empty():
            self._log("DEBUG", "Nucleotides is empty: skipping ahead")
            position = reads[0].reference_start - 1
            contig = reads[0].reference_name
        return (contig, position + 1)

    def _get_column(self, position, bases, region):
        strand = bases.get_strand(threshold=self._strand_confidence_threshold)
        if self._use_strand_correction:
            bases.filter_by_strand(strand)
        if strand == "-":
            bases.complement()

        if not self._valid_column(
                position,
                bases,
                region):
            return None

        variants = bases.get_variants()
        return [
            position + 1, # 1 indexed
            bases.ref,
            strand,
            len(bases),
            f'{bases.mean_quality():.2f}',
            bases.frequencies(),
            ' '.join(variants) if variants else '-',
            f'{bases.get_max_ratio():.2f}',
            '\t'.join(['-'] * 5)
        ]


    def analyze(self, bam_files, region=None):
        '''
        Compute editting rates for a given SAM file.

        Parameters:
            bam_file (str): Path to the SAM file
            region (dict): Region to process. Must have a "contig" key and optional "start"
            and "end".
        '''
        if region is None:
            region = {}

        samfile = AlignmentManager(
            min_quality=self._min_read_quality,
            min_length=self._min_read_length,
            files=bam_files,
            ignore_truncation=True)

        # Open the iterator
        self._log('INFO', 'Fetching data from bams {} [REGION={}]', bam_files, region)
        read_iter = samfile.fetch_by_position(**region)
        reads = next(read_iter, None)

        contig = None
        position = None
        nucleotides = CompiledReads(
            self._strand,
            self._min_base_position,
            self._max_base_position,
            self._min_base_quality)
        if self.reference:
            nucleotides.add_reference(self.reference)
        total = 0

        while reads is not None or not nucleotides.is_empty():
            contig, position = self._next_position(reads, nucleotides, contig, position)
            if position >= region.get("stop", position + 1):
                break
            self._log('DEBUG', 'Analyzing position {} {}', contig, position)

            # Get all the read(s) starting at position
            if reads is not None and reads[0].reference_start == position:
                self._log('DEBUG', 'Adding {} reads', len(reads))
                total += len(reads)
                nucleotides.add_reads(reads)
                reads = next(read_iter, None)

            # Process edits
            bases = nucleotides.pop(position)
            if position in self._exclude_positions.get(contig, []):
                self._log('DEBUG', 'Listed for exclusion - skipping')
                continue
            if self._target_positions and position not in self._target_positions.get(contig, []):
                self._log('DEBUG', 'Not listed for inclusion - skipping')
                continue
            if bases is None:
                self._log('DEBUG', 'No reads - skipping')
                continue
            column = self._get_column(position, bases, region)
            if column is None:
                self._log('DEBUG', 'Bad column - skipping')
                continue
            self._log('DEBUG', 'Yielding output for {} reads', len(bases))
            yield [contig] + column

        self._log('INFO', '[REGION={}] {} total reads read',
                  region,
                  total)

    def set_logging_level(self, level):
        '''
        Set the class logging level.

        The options are "debug", "verbose", or "none".

        Parameters:
            level (str): logging level
        '''
        if level.lower() == "debug":
            self._log = self._log_all
        elif level.lower() == "verbose":
            self._log = self._log_verbose
        else:
            self._log = REDItools._default(None)

    def set_read_filters(self, min_quality=None, min_length=None):
        '''
        Add additional filters to what reads to process. You may set one or both of the arguments.

        Parameters:
            min_quality (int): The minimum quality of reads for inclusion
            min_length (int): The minimum length of reads for inclusion
        '''
        self._min_read_quality = min_quality
        self._min_read_length = min_length

        if self._min_read_quality > 0:
            self._column_checks.add(self._check_column_quality)
        else:
            self._column_checks.discard(self._check_column_quality)

    def set_base_filters(self, min_position=0, max_position=float("inf"), min_quality=0):
        '''
        Add filters to individual bases within reads. You can set each filter individually, or
        all together.

        Parameters:
            min_position (int): The minimum position a base can be within the read (inclusive)
            maX_position (int): The maximum position a base can be within the read
            min_qualtiy (int): The minimum read quality at the base position
        '''
        self._min_base_position = min_position
        self._max_base_position = max_position
        self._min_base_quality = min_quality

    def set_strand(self, strand):
        '''
        Limit the analysis to focus on a specific strand.

        Parameters:
            strand (int): 0 for unstranded, 1 for positive, 2 for negative.
        '''
        self._strand = strand

    def set_strand_confidence_threshold(self, threshold):
        '''
        When a position has mixed strand content, you may set the threshold of certainty for
        which strand to report in output. If the threshold isn't passed, "*" is used.

        Parameters:
            threshold (float): Threshold ratio for determining strandedness
        '''
        self._strand_confidence_threshold = threshold

    def set_min_depth(self, threshold):
        '''
        Only report positions that have a minimum read depth.

        Parameters:
            threshold (int): Minimum number of reads for a given position
        '''
        self._min_column_length = threshold
        if threshold > 1:
            self._column_checks.add(self._check_column_min_length)
        else:
            self._column_checks.discard(self._check_column_min_length)

    def use_strand_correction(self):
        '''
        Only reports reads/positions that match `strand`
        '''
        self._use_strand_correction = True

    def set_min_edits(self, threshold):
        '''
        Only report positions that have a minimum number of bases different from the reference.

        Parameters:
            threshold (int): Minimum number of edits.
        '''
        self._min_edits = threshold
        if threshold > 0:
            self._column_checks.add(self._check_column_edit_frequency)
        else:
            self._column_checks.discard(self._check_column_edit_frequency)

    def set_min_edits_per_nucleotide(self, threshold):
        '''
        Only report positions where all non-reference nucleotides have at least `threshold` reads.

        Parameters:
            threshold (int): Minimum number of edits for all nucleotides
        '''
        self._min_edits_per_nucleotide = threshold
        if threshold > 0:
            self._column_checks.add(self._check_column_min_edits)
        else:
            self._column_checks.discard(self._check_column_min_edits)

    def add_reference(self, reference_fname):
        '''
        Use a reference fasta file instead of reference from the BAM files.
        '''
        self.reference = FastaFile(reference_fname)

    def __init__(self):
        self._hostname_string = utils.get_hostname_string()
        self._min_column_length = 1
        self._min_edits = 0
        self._min_edits_per_nucleotide = 0

        self._log = REDItools._default(None)

        self._strand = 0
        self._use_strand_correction = False
        self._strand_confidence_threshold = 0.5

        self._min_base_quality = 30
        self._min_base_position = 0
        self._max_base_position = float("inf")

        self._column_checks = set()

        self._min_read_length = 30
        self._min_read_quality = 0

        self._target_positions = False
        self._exclude_positions = {}
        self._splice_positions = []
        self._omopolymeric_positions = []

        self.reference = None

class REDItoolsDNA(REDItools):
    '''
    Analysis system for RNA editing events in DNA.
    '''
    def __init__(self):
        self.get_position_strand = REDItools._default("*")
        self._get_strand = REDItools._default("*")
        REDItools.__init__(self)
    def set_strand(self, strand):
        raise ValueError("Cannot set strand value if DNA is True")
