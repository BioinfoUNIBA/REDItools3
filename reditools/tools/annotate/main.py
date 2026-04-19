import csv
import sys
import traceback

import pysam

from reditools import file_utils
from reditools.rtannotater import RTAnnotater
from reditools.tools.annotate.parse_args import parse_args

_contig = 'Region'

def contig_order_from_bam(bam_fname: str) -> dict[str, int]:
    contigs = {}
    with pysam.AlignmentFile(bam_fname, ignore_truncation=True) as bam:
        for idx, contig in enumerate(bam.references, start=1):
            contigs[contig] = idx
    return contigs

def contig_order_from_fai(fai_fname: str) -> dict[str, int]:
    contigs = {}
    with file_utils.open_stream(fai_fname, 'r') as stream:
        for idx, line in enumerate(stream, start=1):
            contig = line.split('\t')[0]
            contigs[contig] = idx
    return contigs

def contig_order_from_out(out_fname: str) -> dict[str, int]:
    contigs: dict[str, int] = {}
    with file_utils.open_stream(out_fname, 'r') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        last_contig = None
        for row in reader:
            if row[_contig] != last_contig:
                if row[_contig] in contigs:
                    raise ValueError(
                        f'File {out_fname} does not appear to be in sorted '
                        'order.'
                    )
                contigs[row[_contig]] = len(contigs) + 1
                last_contig = row[_contig]
    return contigs

def main() -> None:
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

    rta = RTAnnotater(options.rna_file, options.dna_file, contig_order)
    try:
        rta.annotate(sys.stdout)
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
        sys.exit(1)
