import argparse
from reditools import file_utils
import csv
import sys
import traceback
import pysam
from reditools.rtannotater import RTAnnotater
from .parse_args import parse_args

__all__ = ('main')


def contig_order_from_bam(bam_fname):
    contigs = {}
    with pysam.AlignmentFile(bam_fname, ignore_truncation=True) as bam:
        for idx, contig in enumerate(bam.references, start=1):
            contigs[contig] = idx
    return contigs


def contig_order_from_fai(fai_fname):
    contigs = {}
    with file_utils.open_stream(fai_fname, 'r') as stream:
        for idx, line in enumerate(stream, start=1):
            contig = line.split('\t')[0]
            contigs[contig] = idx
    return contigs


def contig_order_from_out(out_fname):
    contigs = {}
    with file_utils.open_stream(out_fname, 'r') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        last_contig = None
        for row in reader:
            if row['Region'] != last_contig:
                if row['Region'] in contigs:
                    raise ValueError(
                        f'File {out_fname} does not appear to be in sorted '
                        'order.'
                    )
                contigs[row['Region']] = len(contigs) + 1
                last_contig = row['Region']
    return contigs


def parse_options():
    """
    Parse commandline options for REDItools.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(
        prog='reditools annotate',
        description='Annotates RNA REDItools output with DNA output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'rna_file',
        help='The REDItools output from RNA data',
    )
    parser.add_argument(
        'dna_file',
        help='The REDItools output from corresponding DNA data',
    )
    parser.add_argument(
        '-d',
        '--debug',
        help='Report stack trace on crash.',
        action='store_true',
    )
    order_group = parser.add_argument_group(
        title='Contig order options',
        description=(
            'By default, the annotate tool will determine the contig order '
            'by reading through the rna_file once before performing the '
            'annotation. Depending on the size of the rna_file, this could '
            'take a while. To skip this time consuming step, you can provide '
            'the contig order from another source.'
        ),
    )
    order_group.add_argument(
        '-b',
        '--bam',
        help='BAM file to get contig order from.',
    )
    order_group.add_argument(
        '-f',
        '--fai',
        help='FASTA Index file to get contig order from.',
    )
    options = parser.parse_args()

    if options.bam and options.fai:
        parser.error(
            message='Options -b/--bam and -f/--fai are mutually exclusive.',
        )

    return options


def main():
    options = parse_args()
    try:
        if options.fai:
            order_fname = options.fai
            contig_order = contig_order_from_fai(options.fai)
        elif options.bam:
            order_fname = options.bam
            contig_order = contig_order_from_bam(options.bam)
        else:
            order_fname = options.rna_file
            contig_order = contig_order_from_out(options.rna_file)
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(
            '[ERROR] There was an error getting the contig order from '
            f'{order_fname}: ({type(exc)}) {exc}\n'
        )
        sys.exit(1)

    try:
        x = RTAnnotater(options.rna_file, options.dna_file, contig_order)
        x.annotate(sys.stdout)
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
        sys.exit(1)
