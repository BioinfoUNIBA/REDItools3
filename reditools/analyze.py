"""Commandline tool for REDItools."""

import argparse
import csv
import sys
import traceback
from multiprocessing import Process, Queue
from queue import Empty as EmptyQueueException
from tempfile import NamedTemporaryFile

from reditools import file_utils, reditools, utils
from reditools.alignment_manager import AlignmentManager
from reditools.logger import Logger
from reditools.region import Region

_contig = 'contig'
_start = 'start'
_stop = 'stop'

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


def setup_alignment_manager(options):
    """
    Create an AlignmentManager for REDItools.

    Parameters:
        options (namespace): Commandline arguments

    Returns:
        AlignmentManager
    """
    sam_manager = AlignmentManager(
        ignore_truncation=True,
    )
    sam_manager.min_quality = options.min_read_quality
    sam_manager.min_length = options.min_read_length
    for sam in options.file:
        sam_manager.add_file(
            sam,
            options.exclude_reads,
        )
    return sam_manager


def setup_rtools(options):  # noqa:WPS213,WPS231
    """
    Create a REDItools object.

    Parameters:
        options (namespace): Commandline arguments from argparse

    Returns:
        A configured REDItools object
    """
    if options.dna:
        rtools = reditools.REDItoolsDNA()
    else:
        rtools = reditools.REDItools()

    if options.debug:
        rtools.log_level = Logger.debug_level
    elif options.verbose:
        rtools.log_level = Logger.info_level

    if options.load_omopolymeric_file:
        regions = file_utils.load_omopolymeric_regions(
            options.load_omopolymeric_file,
        )
        rtools.exclude(regions)

    if options.create_omopolymeric_file:
        rtools.create_omopolymeric_positions(
            options.create_omopolymeric_file,
            options.omopolymeric_span,
        )

    if options.splicing_file:
        rtools.load_splicing_file(
            options.splicing_file,
            options.splicing_span,
        )

    if options.bed_file:
        rtools.load_target_positions(options.bed_file)
    if options.exclude_regions:
        for fname in options.exclude_regions:
            regions = file_utils.read_bed_file(fname)
            rtools.exclude(regions)
    if options.reference:
        rtools.add_reference(options.reference)

    rtools.min_base_position = options.min_base_position
    rtools.max_base_position = options.max_base_position
    rtools.min_base_quality = options.min_base_quality

    rtools.min_column_length = options.min_column_length
    rtools.min_edits = options.min_edits
    rtools.min_edits_per_nucleotide = options.min_edits_per_nucleotide
    rtools.strand = options.strand

    rtools.strand_confidence_threshold = options.strand_confidence_threshold

    if options.strand_correction:
        rtools.use_strand_correction()
    if options.exclude_multis:
        rtools.only_one_alt()

    return rtools


def region_args(bam_fname, region, window):
    """
    Split a region into segments for paralllel processing.

    Parameters:
        bam_fname (str): BAM file to collect contig info from
        region (Region): Genomic region to split
        window (int): How large the sub regions should be.

    Returns:
        (list): Sub regions
    """
    if region is not None:
        if window:
            return region.split(window)
        return [region]

    args = []
    for contig, size in utils.get_contigs(bam_fname):
        region = Region(contig=contig, start=1, stop=size+1)
        if window:
            args.extend(region.split(window))
        else:
            args.append(region)
    return args


def write_results(rtools, sam_manager, file_name, region, output_format):
    """
    Write the results from a REDItools analysis to a temporary file.

    Parameters:
        rtools (REDItools): REDItools instance
        sam_manager (AlignmentManager): Source of reads
        file_name (string): Input file name for analysis
        region: Region to analyze
        output_format (dict): keyword arguments for csv.writer constructor.

    Returns:
        string: Name of the temporary file.
    """
    with NamedTemporaryFile(mode='w', delete=False) as stream:
        writer = csv.writer(stream, **output_format)
        for rt_result in rtools.analyze(sam_manager, region):
            variants = rt_result.variants
            writer.writerow([
                rt_result.contig,
                rt_result.position,
                rt_result.reference,
                rt_result.strand,
                rt_result.depth,
                f'{rt_result.mean_quality:.2f}',
                rt_result.per_base_depth,
                ' '.join(sorted(variants)) if variants else '-',
                f'{rt_result.edit_ratio:.2f}',
                '\t'.join(['-' for _ in range(5)]),
            ])
        return stream.name


def run(options, in_queue, out_queue):
    """
    Analyze a genomic segment using REDItools.

    Parameters:
        options (namesapce): Configuration options from argparse for REDItools
        in_queue (Queue): Queue of input arguments for analysis
        out_queue (Queue): Queue to store paths to analysis results

    Returns:
        bool: True if the in_queue is empty
    """
    try:
        rtools = setup_rtools(options)
        while True:
            args = in_queue.get()
            if args is None:
                return True
            sam_manager = setup_alignment_manager(options)
            idx, region = args
            file_name = write_results(
                rtools,
                sam_manager,
                options.file,
                region,
                options.output_format,
            )
            out_queue.put((idx, file_name))
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] {exc}\n')


def parse_options():  # noqa:WPS213
    """
    Parse commandline options for REDItools.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument(
        'file',
        nargs='+',
        help='The bam file to be analyzed',
    )
    parser.add_argument(
        '-r',
        '--reference',
        help='The reference FASTA file',
    )
    parser.add_argument(
        '-o',
        '--output-file',
        help='The output statistics file',
    )
    parser.add_argument(
        '-s',
        '--strand',
        choices=(0, 1, 2),
        type=int,
        default=0,
        help='Strand: this can be 0 (unstranded),' +
        '1 (secondstrand oriented) or ' +
        '2 (firststrand oriented)',
    )
    parser.add_argument(
        '-a',
        '--append-file',
        action='store_true',
        help='Appends results to file (and creates if not existing)',
    )
    parser.add_argument(
        '-g',
        '--region',
        help='The self.region of the bam file to be analyzed',
    )
    parser.add_argument(
        '-m',
        '--load-omopolymeric-file',
        help='The file containing the omopolymeric positions',
    )
    parser.add_argument(
        '-c',
        '--create-omopolymeric-file',
        default=False,
        help='Path to write omopolymeric positions to',
        action='store_true',
    )
    parser.add_argument(
        '-os',
        '--omopolymeric-span',
        type=int,
        default=5,
        help='The omopolymeric span',
    )
    parser.add_argument(
        '-sf',
        '--splicing-file',
        help='The file containing the splicing sites positions',
    )
    parser.add_argument(
        '-ss',
        '--splicing-span',
        type=int,
        default=4,
        help='The splicing span',
    )
    parser.add_argument(
        '-mrl',
        '--min-read-length',
        type=int,
        default=30,  # noqa:WPS432
        help='Reads whose length is below this value will be discarded.',
    )
    parser.add_argument(
        '-q',
        '--min-read-quality',
        type=int,
        default=20,  # noqa:WPS432
        help='Reads with mapping quality below this value will be discarded.',
    )
    parser.add_argument(
        '-bq',
        '--min-base-quality',
        type=int,
        default=30,  # noqa:WPS432
        help='Base quality below this value will not be included in ' +
        'the analysis.',
    )
    parser.add_argument(
        '-mbp',
        '--min-base-position',
        type=int,
        default=0,
        help='Bases which reside in a previous position (in the read)' +
        'will not be included in the analysis.',
    )
    parser.add_argument(
        '-Mbp',
        '--max-base-position',
        type=int,
        default=0,
        help='Bases which reside in a further position (in the read)' +
        'will not be included in the analysis.',
    )
    parser.add_argument(
        '-l',
        '--min-column-length',
        type=int,
        default=1,
        help='Positions whose columns have length below this value will' +
        'not be included in the analysis.',
    )
    parser.add_argument(
        '-e',
        '--exclude-multis',
        default=False,
        help='Do not report any position with more than one alternate base.',
        action='store_true',
    )
    parser.add_argument(
        '-men',
        '--min-edits-per-nucleotide',
        type=int,
        default=0,
        help='Positions whose columns have bases with less than' +
        'min-edits-per-base edits will not be included in the analysis.',
    )
    parser.add_argument(
        '-me',
        '--min-edits',
        type=int,
        default=0,  # noqa:WPS432
        help='The minimum number of editing events (per position). ' +
        'Positions whose columns have bases with less than ' +
        '"min-edits-per-base edits" will not be included in the ' +
        'analysis.',
    )
    parser.add_argument(
        '-Men',
        '--max-editing-nucleotides',
        type=int,
        default=100,  # noqa:WPS432
        help='The maximum number of editing nucleotides, from 0 to 4 ' +
        '(per position). Positions whose columns have more than ' +
        '"max-editing-nucleotides" will not be included in the analysis.',
    )
    parser.add_argument(
        '-T',
        '--strand-confidence-threshold',
        type=float,
        default=0.7,  # noqa:WPS432
        help='Only report the strandedness if at least this proportion of ' +
        'reads are of a given strand',
    )
    parser.add_argument(
        '-C',
        '--strand-correction',
        default=False,
        help='Strand correction. Once the strand has been inferred, ' +
        'only bases according to this strand will be selected.',
        action='store_true',
    )
    parser.add_argument(
        '-V',
        '--verbose',
        default=False,
        help='Verbose information in stderr',
        action='store_true',
    )
    parser.add_argument(
        '-N',
        '--dna',
        default=False,
        help='Run REDItools 2.0 on DNA-Seq data',
        action='store_true',
    )
    parser.add_argument(
        '-B',
        '--bed_file',
        help='Path of BED file containing target self.regions',
    )
    parser.add_argument(
        '-t',
        '--threads',
        help='Number of threads to run',
        type=int,
        default=1,
    )
    parser.add_argument(
        '-w',
        '--window',
        help='How many bp should be processed by each thread at a time. ' +
        'Defaults to full contig.',
        type=int,
        default=0,
    )
    parser.add_argument(
        '-k',
        '--exclude_regions',
        nargs='+',
        help='Path of BED file containing regions to exclude from analysis',
    )
    parser.add_argument(
        '-E',
        '--exclude_reads',
        help='Path to a text file listing read names to exclude from analysis',
    )
    parser.add_argument(
        '-d',
        '--debug',
        default=False,
        help='REDItools is run in DEBUG mode.',
        action='store_true',
    )

    return parser.parse_args()


def check_dead(processes):
    """
    Look through processes to determine if any have died unexpectedly.

    If any process has an exit code of 1, this method will terminate all other
    processes and then exit with code 1.

    Parameters:
        processes (list): Processes to check
    """
    for proc in processes:
        if proc.exitcode == 1:
            for to_kill in processes:
                to_kill.kill()
            sys.stderr.write('[ERROR] Killing job\n')
            sys.exit(1)


def main():
    """Perform RNA editing analysis."""
    options = parse_options()
    options.output_format = {'delimiter': '\t', 'lineterminator': '\n'}
    options.encoding = 'utf-8'
    if options.exclude_reads:
        options.exclude_reads = file_utils.load_text_file(
            options.exclude_reads,
        )

    # Put analysis chunks into queue
    regions = region_args(
        options.file[0],
        Region(string=options.region) if options.region else None,
        window=options.window,
    )

    in_queue = Queue()
    for args in enumerate(regions):
        in_queue.put(args)
    for _ in range(options.threads):
        in_queue.put(None)

    # Start parallel jobs
    out_queue = Queue()
    processes = [
        Process(
            target=run,
            args=(options, in_queue, out_queue),
        ) for _ in range(options.threads)
    ]
    concat_output(
        options,
        monitor(processes, out_queue, in_queue.qsize()),
    )


def monitor(processes, out_queue, chunks):
    """
    Monitor parallel REDItools jobs.

    Parameters:
        processes (list): Threads
        out_queue (Queue): Output of threads
        chunks (int): Number of chunks for analysis

    Returns:
        list: Temporary files containing the output of each chunk.
    """
    tfs = [None for _ in range(chunks - len(processes))]

    for prc in processes:
        prc.start()

    while None in tfs:
        try:
            idx, fname = out_queue.get(block=False, timeout=1)
            tfs[idx] = fname
        except EmptyQueueException:
            check_dead(processes)
    return tfs


def concat_output(options, tfs):
    """
    Write the output of a REDItools analysis.

    Parameters:
        options (namespace): Commandline options for file formatting.
        tfs (list): Temporary files containing REDItools results
    """
    # Setup final output file
    if options.output_file:
        mode = 'a' if options.append_file else 'w'
        stream = file_utils.open_stream(
            options.output_file,
            mode,
            encoding=options.encoding,
        )
    else:
        stream = sys.stdout

    with stream:
        writer = csv.writer(stream, **options.output_format)
        if not options.append_file:
            writer.writerow(fieldnames)
        file_utils.concat(stream, *tfs, encoding=options.encoding)
