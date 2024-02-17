#!/usr/bin/env python
'''
Commandline tool for REDItools
'''
from __future__ import division

import argparse
import csv
import re
import os
import sys
import traceback

from queue import Empty as EmptyQueueException
from tempfile import NamedTemporaryFile
from multiprocessing import Process
from multiprocessing import Queue
from pysam.libcalignmentfile import AlignmentFile

from reditools import reditools


def setup(options):
    '''Creates a REDItools object.

    Parameters:
        options (namespace): Commandline arguments from argparse

    Returns:
        A configured REDItools object
    '''

    if options.dna:
        rtools = reditools.REDItoolsDNA()
    else:
        rtools = reditools.REDItools()

    if options.debug:
        rtools.set_logging_level("debug")
    elif options.verbose:
        rtools.set_logging_level("verbose")

    if options.load_omopolymeric_file:
        rtools.load_omopolymeric_positions(options.load_omopolymeric_file)

    if options.create_omopolymeric_file:
        rtools.create_omopolymeric_positions(
            options.create_omopolymeric_file,
            options.omopolymeric_span)

    if options.splicing_file:
        rtools.load_splicing_file(
            options.splicing_file,
            options.splicing_span)

    if options.bed_file:
        rtools.load_target_positions(options.bed_file)
    if options.exclude_regions:
        rtools.load_exclude_positions(options.exclude_regions)

    if options.reference:
        rtools.add_reference(options.reference)

    rtools.set_read_filters(
        min_quality=options.min_read_quality,
        min_length=options.min_read_length)

    rtools.set_base_filters(
        min_position=options.min_base_position,
        max_position=options.max_base_position,
        min_quality=options.min_base_quality)

    rtools.set_min_depth(options.min_column_length)
    rtools.set_min_edits(options.min_edits)
    rtools.set_min_edits_per_nucleotide(options.min_edits_per_nucleotide)
    rtools.set_strand(options.strand)

    rtools.set_strand_confidence_threshold(options.strand_confidence_threshold)

    if options.strand_correction:
        rtools.use_strand_correction()

    return rtools


def contig_window_args(contig, start, window, end, i=0):
    '''Produces regional segments based on a window size.

    Parameters:
        contig (string): Contig or chromsome name
        start (int): Region start position
        window (int): Window size in bp
        end (int): Region end
        i (int): Region order for recombining parallel processing output

    Returns:
        Generator of tuples (i, region)
    '''
    while start + window < end:
        yield (i, {"contig": contig, "start": start, "stop": start + window})
        start += window
        i += 1
    yield (i, {"contig": contig, "start": start, "stop": end})

def window_args(contigs, sizes, window):
    '''Produces region segments for a genome.

    Parameters:
        contigs (iterable): Contigs or chromsome names
        sizes (iterable): The sizes of each contig or chromsome (must be same length and
            order as `contigs`)
        window (int): Window size in bp

    Returns:
        Generator of tuples (i, region)
    '''
    i = 0
    for contig, size in zip(contigs, sizes):
        for arg in contig_window_args(contig, 0, window, size, i):
            yield arg
            i += 1

def region_args(region, window):
    '''Produces regional segments based on a window size.

    Parameters:
        region (dict): Genomic region with keys "contig", "start", and "stop"
        window (int): Window size in bp
        i (int): Region order for recombining parallel processing output

    Returns:
        Generator of tuples (i, region)
    '''

    return contig_window_args(region["contig"], region["start"], window, region["stop"])

def get_args(options):
    '''Produces arguments to `run` for parallel processing.

    Parameters:
        options (namespace): Command line options from parseargs

    Returns:
        Generator of arguments for `run`
    '''
    region = parse_region(options.region) if options.region else {}
    contigs, sizes = get_contigs(options.file[0])

    # Put analysis chunks into queue
    if options.window:
        if region:
            size = sizes[contigs.index(region["contig"])]
            return region_args({
                "contig": region["contig"],
                "start": region.get("start", 0),
                "stop": region.get("stop", size)}, options.window)
        return window_args(contigs, sizes, options.window)
    if region:
        return [(0, region)]
    return ((i, {"contig": contigs[i]}) for i in range(len(contigs)))

def run(options, in_queue, out_queue):
    '''Analyzes a genomic segment using REDItools.

    Parameters:
        options (namesapce): Configuration options from argparse for REDItools
        in_queue (Queue): Queue of input arguments for analysis
        out_queue (Queue): Queue to store paths to analysis results
    '''
    try:
        rtools = setup(options)
        while True:
            args = in_queue.get()
            if args is None:
                return True
            i, region = args
            with NamedTemporaryFile(mode="w", delete=False) as stream:
                writer = csv.writer(stream, **options.output_format)
                for result in rtools.analyze(options.file, region):
                    # Transformations for REDItools2 comparison
                    result = result[:-1] + result[-1].split()
                    result[3] = {"-": 0, "+": 1, "*": 2}[result[3]]
                    writer.writerow(result)
                out_queue.put((i, stream.name))
    except Exception as exc: # pylint: disable=broad-except
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f"[ERROR] {exc}\n")
        sys.exit(1)

def to_int(string):
    '''Converts a potentially formatted string to an int

    Parameters:
        string (str): A string representation of an integer

    Returns:
        The integer values of the string.
    '''
    return int(re.sub(r"[\s,]", "", string))

def parse_region(region_str):
    '''
    Parses a region string into chromosome, start, and end.

    Parameters:
        region_str (str): In the format of 'chr#' or 'chr#:start-end'

    Returns:
        A dict representation of the region with keys "contig", "start", and "stop".
    '''
    if region_str is None:
        return None
    region = re.split("[:-]", region_str)
    if region and len(region) in (1, 3):
        if len(region) == 1:
            return {"contig": region[0]}
        if len(region) == 3:
            start = to_int(region[1])
            stop = to_int(region[2])
            if start >= stop:
                raise ValueError("Please provide a region of the form chrom:start-end " +
                                f"(with end > start). Region provided: {region}")
            return {"contig": region[0], "start": start, "stop": stop}
    raise ValueError("Please provide a region of the form chrom:start-end " +
                    f"(with end > start). Region provided: {region}")

def parse_options():
    '''
    Parses commandline options for REDItools.
    '''
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('file',
         nargs='+',
		 help='The bam file to be analyzed')
    parser.add_argument('-r',
         '--reference',
		 help='The reference FASTA file')
    parser.add_argument('-o',
		 '--output-file',
		 help='The output statistics file')
    parser.add_argument('-s',
		 '--strand',
         choices=(0, 1, 2),
		 type=int,
		 default=0,
		 help='Strand: this can be 0 (unstranded),' +
            '1 (secondstrand oriented) or ' +
            '2 (firststrand oriented)')
    parser.add_argument('-a',
		 '--append-file',
		 action='store_true',
		 help='Appends results to file (and creates if not existing)')
    parser.add_argument('-g',
		 '--region',
		 help='The self.region of the bam file to be analyzed')
    parser.add_argument('-m',
		 '--load-omopolymeric-file',
		 help='The file containing the omopolymeric positions')
    parser.add_argument('-c',
		 '--create-omopolymeric-file',
		 default=False,
		 help='Path to write omopolymeric positions to',
		 action='store_true')
    parser.add_argument('-os',
		 '--omopolymeric-span',
		 type=int,
		 default=5,
		 help='The omopolymeric span')
    parser.add_argument('-sf',
		 '--splicing-file',
		 help='The file containing the splicing sites positions')
    parser.add_argument('-ss',
		 '--splicing-span',
		 type=int,
		 default=4,
		 help='The splicing span')
    parser.add_argument('-mrl',
		 '--min-read-length',
		 type=int,
		 default=30,
		 help='Reads whose length is below this value will be discarded.')
    parser.add_argument('-q',
		 '--min-read-quality',
		 type=int,
		 default=20,
		 help='Reads whose mapping quality is below this value will be discarded.')
    parser.add_argument('-bq',
		 '--min-base-quality',
		 type=int,
		 default=30,
		 help='Bases whose quality is below this value will not be included in the analysis.')
    parser.add_argument('-mbp',
		 '--min-base-position',
		 type=int,
		 default=0,
		 help='Bases which reside in a previous position (in the read)' +
            'will not be included in the analysis.')
    parser.add_argument('-Mbp',
		 '--max-base-position',
		 type=int,
		 default=0,
		 help='Bases which reside in a further position (in the read)' +
            'will not be included in the analysis.')
    parser.add_argument('-l',
		 '--min-column-length',
		 type=int,
		 default=1,
		 help='Positions whose columns have length below this value will' +
            'not be included in the analysis.')
    parser.add_argument('-men',
		 '--min-edits-per-nucleotide',
		 type=int,
		 default=0,
		 help='Positions whose columns have bases with less than' +
            'min-edits-per-base edits will not be included in the analysis.')
    parser.add_argument('-me',
		 '--min-edits',
		 type=int,
		 default=0,
		 help='The minimum number of editing events (per position). ' +
            'Positions whose columns have bases with less than ' +
            '"min-edits-per-base edits" will not be included in the analysis.')
    parser.add_argument('-Men',
		 '--max-editing-nucleotides',
		 type=int,
		 default=100,
		 help='The maximum number of editing nucleotides, from 0 to 4 (per position). ' +
            'Positions whose columns have more than "max-editing-nucleotides" ' +
            'will not be included in the analysis.')
    parser.add_argument('-d',
		 '--debug',
		 default=False,
		 help='REDItools is run in DEBUG mode.',
		 action='store_true')
    parser.add_argument('-T',
		 '--strand-confidence-threshold',
         type=float,
		 default=0.7,
		 help='Only report the strandedness if at least this proportion of reads are of a given strand')
    parser.add_argument('-C',
		 '--strand-correction',
		 default=False,
		 help='Strand correction. Once the strand has been inferred, ' +
            'only bases according to this strand will be selected.',
		 action='store_true')
    parser.add_argument('-V',
		 '--verbose',
		 default=False,
		 help='Verbose information in stderr',
		 action='store_true')
    parser.add_argument('-N',
		 '--dna',
		 default=False,
		 help='Run REDItools 2.0 on DNA-Seq data',
		 action='store_true')
    parser.add_argument('-B',
		 '--bed_file',
		 help='Path of BED file containing target self.regions')
    parser.add_argument('-t',
        '--threads',
        help='Number of threads to run',
        type=int,
        default=1)
    parser.add_argument('-w',
        '--window',
        help='How many bp should be processed by each thread at a time. Defaults to full contig.',
        type=int,
        default=0)
    parser.add_argument('-k',
        '--exclude_regions',
        help='Path of BED file containing regions to exclude from analysis')

    args = parser.parse_args()

    return args

def get_contigs(sam_path):
    '''Retrieves contig or chromsome data from an alignment file.

    Parameters:
        sam_path (string): Path to an alignment file.

    Returns:
        tuple of lists containing the reference names and reference lengths in corresponding order
    '''
    with AlignmentFile(sam_path, ignore_truncation=True) as sam:
        contigs = list(sam.references)
        sizes = list(sam.lengths)
        indices = sorted(list(range(len(contigs))), key=lambda i: contigs[i])
        return ([contigs[i] for i in indices], [sizes[i] for i in indices])

def concat(output, *fnames, clean_up=True, encoding="utf-8"):
    '''Combines one or more files into another file.

    Parameters:
        output (file): A file like object for writing
        *fnames (string): Paths to files for concatenation
        clean_up (bool): If True, deletes the files in *fnames after concatenation
        encoding (string): File encoding
    '''
    for fname in fnames:
        with open(fname, "r", encoding=encoding) as stream:
            for line in stream:
                output.write(line)
        if clean_up:
            os.remove(fname)

def check_dead(processes):
    '''
    Looks through a list of processes to determine if any have died unexpectedly.
    If any process has an exit code of 1, this method will terminate all other
    processes and then exit with code 1.

    Parameters:
        processes (list) - Processes to check
    '''
    for proc in processes:
        if proc.exitcode == 1:
            for to_kill in processes:
                to_kill.kill()
            sys.stderr.write("[ERROR] Killing job\n")
            sys.exit(1)

def main():
    '''Performs RNA editing analysis'''
    options = parse_options()

    options.output_format = {"delimiter": "\t", "lineterminator": "\n"}
    options.encoding = "utf-8"

    # Put analysis chunks into queue
    in_queue = Queue()
    num_args = 0
    for arg in get_args(options):
        in_queue.put(arg)
        num_args += 1
    for i in range(options.threads):
        in_queue.put(None)

    # Start parallel jobs
    processes = []
    out_queue = Queue()
    for i in range(options.threads):
        processes.append(Process(target=run, args=(options, in_queue, out_queue)))
        processes[-1].start()

    tfs = [None] * num_args
    i = 0
    while i < num_args:
        try:
            j, fname = out_queue.get(block=False, timeout=1)
            tfs[j] = fname
            i += 1
        except EmptyQueueException:
            check_dead(processes)


    # Setup final output file
    if options.output_file:
        mode = "a" if options.append_file else "w"
        stream = open(options.output_file, mode, encoding=options.encoding)
    else:
        stream = sys.stdout

    with stream:
        writer = csv.writer(stream, **options.output_format)
        if not options.append_file:
            writer.writerow(reditools.REDItools.fieldnames)
        concat(stream, *tfs, encoding=options.encoding)

if __name__ == '__main__':
    main()
