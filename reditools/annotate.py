import argparse
from reditools import file_utils
import csv
import sys


class RTAnnotater:
    def __init__(self, rna_file, dna_file):
        self.rna_file = rna_file
        self.dna_file = dna_file
        self.contig_order = self._load_contig_order()

    def _load_contig_order(self):
        contigs = {}
        idx = 1
        with file_utils.open_stream(self.rna_file, 'r') as stream:
            reader = csv.reader(stream, delimiter='\t')
            last_contig = next(reader)[0]
            contigs[last_contig] = 0
            for row in reader:
                if row[0] == last_contig:
                    continue
                contigs[row[0]] = idx
                idx += 1
                last_contig = row[0]
        return contigs

    def _cmp_position(self, rna_contig, rna_pos, dna_contig, dna_pos):
        rna_contig = self.contig_order[rna_contig]
        dna_contig = self.contig_order.get(dna_contig, len(self.contig_order))
        if rna_contig < dna_contig:
            return -1
        if rna_contig > dna_contig:
            return 1
        rna_pos = int(rna_pos)
        dna_pos = int(dna_pos)
        if rna_pos < dna_pos:
            return -1
        if rna_pos > dna_pos:
            return 1
        return 0

    def _annotate_row(self, rna_row, dna_row):
        rna_row['gCoverage-q30'] = dna_row['Coverage-q30']
        rna_row['gMeanQ'] = dna_row['MeanQ']
        rna_row['gBaseCount[A,C,G,T]'] = dna_row['BaseCount[A,C,G,T]']
        rna_row['gAllSubs'] = dna_row['AllSubs']
        rna_row['gFrequency'] = dna_row['Frequency']
        return rna_row

    def _compare_files(self):
        with file_utils.open_stream(self.rna_file, 'r') as rna_stream, \
                file_utils.open_stream(self.dna_file, 'r') as dna_stream:
            rna_reader = csv.DictReader(rna_stream, delimiter='\t')
            dna_reader = csv.DictReader(dna_stream, delimiter='\t')

            rna_entry = next(rna_reader, None)
            dna_entry = next(dna_reader, None)

            while rna_entry is not None:
                if dna_entry is None:
                    yield rna_entry
                    rna_entry = next(rna_reader, None)
                    continue
                cmp = self._cmp_position(
                    rna_entry['Region'],
                    rna_entry['Position'],
                    dna_entry['Region'],
                    dna_entry['Position'])
                if cmp == 0:
                    yield self._annotate_row(rna_entry, dna_entry)
                    rna_entry = next(rna_reader, None)
                    dna_entry = next(dna_reader, None)
                elif cmp > 0:
                    dna_entry = next(dna_reader, None)
                else:
                    yield rna_entry
                    rna_entry = next(rna_reader, None)

    def annotate(self, stream):
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=[
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
            'gFrequency'])
        writer.writeheader()
        writer.writerows(self._compare_files())


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
    return parser.parse_args()


def main():
    options = parse_options()
    x = RTAnnotater(options.rna_file, options.dna_file)
    x.annotate(sys.stdout)
