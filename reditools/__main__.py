"""Commandline tool for REDItools."""

import argparse
import csv
import re
import sys
import traceback
from multiprocessing import Process, Queue
from queue import Empty as EmptyQueueException
from tempfile import NamedTemporaryFile

from reditools import reditools, utils
from reditools.logger import Logger

_contig = 'contig'
_start = 'start'
_stop = 'stop'


def setup(options):  # noqa:WPS213
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
        rtools.load_omopolymeric_positions(options.load_omopolymeric_file)

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
        rtools.load_exclude_positions(options.exclude_regions)

    if options.reference:
        rtools.add_reference(options.reference)

    rtools.min_quality = options.min_read_quality
    rtools.min_length = options.min_read_length

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

    return rtools


def contig_window_args(contig, start, window, end, idx=0):
    """
    Produce regional segments based on a window size.

    Parameters:
        contig (string): Contig or chromsome name
        start (int): Region start position
        window (int): Window size in bp
        end (int): Region end
        idx (int): Region order for recombining parallel processing output

    Yields:
        tuples (idx, region)
    """
    while start + window < end:
        yield (idx, {_contig: contig, _start: start, _stop: start + window})
        start += window
        idx += 1
    yield (idx, {_contig: contig, _start: start, _stop: end})


def window_args(contigs, sizes, window):
    """
    Produce region segments.

    Parameters:
        contigs (iterable): Contigs/chromsome names
        sizes (iterable): Contig sizes
        window (int): Window size in bp

    Yields:
        tuples (i, region)
    """
    idx = 0
    for contig, size in zip(contigs, sizes):
        for arg in contig_window_args(contig, 0, window, size, idx):
            yield arg
            idx += 1


def region_args(region, window):
    """
    Produce regional segments based on a window size.

    Parameters:
        region (dict): Genomic region with keys "contig", "start", and "stop"
        window (int): Window size in bp

    Returns:
        Generator of tuples (i, region)
    """
    return contig_window_args(
        region[_contig],
        region[_start],
        window,
        region[_stop],
    )


def get_args(options):
    """
    Produce arguments to `run` for parallel processing.

    Parameters:
        options (namespace): Command line options from parseargs

    Returns:
        Generator of arguments for `run`
    """
    region = parse_region(options.region) if options.region else {}
    contigs, sizes = utils.get_contigs(options.file[0])

    # Put analysis chunks into queue
    if options.window:
        if region:
            size = sizes[contigs.index(region[_contig])]
            return region_args(
                {
                    _contig: region[_contig],
                    _start: region.get(_start, 0),
                    _stop: region.get(_stop, size),
                },
                options.window,
            )
        return window_args(contigs, sizes, options.window)
    if region:
        return [(0, region)]
    c_range = range(len(contigs))
    return ((idx, {_contig: contigs[idx]}) for idx in c_range)


def write_results(rtools, file_name, region, output_format):
    """
    Write the results from a REDItools analysis to a temporary file.

    Parameters:
        rtools (REDItools): REDItools instance
        file_name (string): Input file name for analysis
        region: Region to analyze
        output_format (dict): keyword arguments for csv.writer constructor.

    Returns:
        string: Name of the temporary file.
    """
    with NamedTemporaryFile(mode='w', delete=False) as stream:
        writer = csv.writer(stream, **output_format)
        for rtools_output in rtools.analyze(file_name, region):
            writer.writerow(rtools_output)
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
        rtools = setup(options)
        while True:
            args = in_queue.get()
            if args is None:
                return True
            idx, region = args
            file_name = write_results(
                rtools,
                options.file,
                region,
                options.output_format,
            )
            out_queue.put((idx, file_name))
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] {exc}\n')


def parse_region(region_str):
    """
    Parse a region string into chromosome, start, and end.

    Parameters:
        region_str (str): In the format of 'chr#' or 'chr#:start-end'

    Returns:
        A dict of the region with keys "contig", "start", and "stop".

    Raises:
        ValueError: On improper string format
    """
    if region_str is None:
        return None
    region = re.split('[:-]', region_str)
    if not region:
        return None
    if len(region) == 1:
        return {_contig: region[0]}
    if len(region) == 3:
        start = utils.to_int(region[1])
        stop = utils.to_int(region[2])
        if start >= stop:
            raise ValueError(
                'Please provide a region of the form chrom:' +
                f'start-end (with end > start). Region provided: {region}',
            )
        return {_contig: region[0], _start: start, _stop: stop}
    raise ValueError(
        'Please provide a region of the form chrom:start-end ' +
        f'(with end > start). Region provided: {region}',
    )


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
        '-d',
        '--debug',
        default=False,
        help='REDItools is run in DEBUG mode.',
        action='store_true',
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
        help='Path of BED file containing regions to exclude from analysis',
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

    # Put analysis chunks into queue
    in_queue = Queue()
    for arg in get_args(options):
        in_queue.put(arg)
    for _ in range(options.threads):
        in_queue.put(None)

    # Start parallel jobs
    processes = []
    out_queue = Queue()
    for _ in range(options.threads):
        processes.append(
            Process(
                target=run,
                args=(options, in_queue, out_queue),
            ),
        )
    tfs = monitor(processes, out_queue, in_queue.qsize())
    concat_output(options, tfs)


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
        stream = open(  # noqa:WPS515
            options.output_file,
            mode,
            encoding=options.encoding,
        )
    else:
        stream = sys.stdout

    with stream:
        writer = csv.writer(stream, **options.output_format)
        if not options.append_file:
            writer.writerow(reditools.REDItools.fieldnames)
        utils.concat(stream, *tfs, encoding=options.encoding)


if __name__ == '__main__':
    main()
