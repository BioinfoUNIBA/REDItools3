import argparse
from reditools import file_utils
import csv
import sys
import traceback


class RTAnnotater:
    def __init__(self, rna_file, dna_file):
        self.rna_file = rna_file
        self.dna_file = dna_file
        self.contig_order = self._load_contig_order(rna_file)

    def _load_contig_order(self, fname):
        contigs = {}
        with file_utils.open_stream(fname, 'r') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            last_contig = None
            for row in reader:
                if row['Region'] != last_contig:
                    if row['Region'] in contigs:
                        raise ValueError(
                            f'File {fname} does not appear to be in sorted '
                            'order.'
                        )
                    contigs[row['Region']] = len(contigs) + 1
                    last_contig = row['Region']
        return contigs

    def _cmp_position(self, rna_entry, dna_entry):
        if dna_entry is None:
            return -1
        rna_contig_idx = self.contig_order[rna_entry['Region']]
        # If the DNA contig is not in the RNA file, assume its position is
        # earlier than the current RNA contig to induce fast-forwarding.
        dna_contig_idx = self.contig_order.get(
            dna_entry['Region'],
            0
        )
        if rna_contig_idx == dna_contig_idx:
            return int(rna_entry['Position']) - int(dna_entry['Position'])
        return rna_contig_idx - dna_contig_idx

    def _annotate_row(self, rna_row, dna_row):
        rna_row['gCoverage'] = dna_row['Coverage']
        rna_row['gMeanQ'] = dna_row['MeanQ']
        rna_row['gBaseCount[A,C,G,T]'] = dna_row['BaseCount[A,C,G,T]']
        rna_row['gAllSubs'] = dna_row['AllSubs']
        rna_row['gFrequency'] = dna_row['Frequency']
        return rna_row

    def _legacy_translate(self, row, old_key, new_key):
        row[new_key] = row.pop(old_key, row.get(new_key))
        return row

    def _compare_files(self):
        with file_utils.open_stream(self.rna_file, 'r') as rna_stream, \
                file_utils.open_stream(self.dna_file, 'r') as dna_stream:
            rna_reader = csv.DictReader(rna_stream, delimiter='\t')
            dna_reader = csv.DictReader(dna_stream, delimiter='\t')

            dna_entry = next(dna_reader, None)

            for rna_entry in rna_reader:
                self._legacy_translate(rna_entry, 'Coverage-q30', 'Coverage')
                self._legacy_translate(rna_entry, 'gCoverage-q30', 'gCoverage')

                pos_comp = self._cmp_position(rna_entry, dna_entry)
                while pos_comp > 0:
                    dna_entry = next(dna_reader, None)
                    pos_comp = self._cmp_position(rna_entry, dna_entry)
                if pos_comp == 0:
                    self._legacy_translate(
                        dna_entry,
                        'Coverage-q30',
                        'Coverage',
                    )
                    yield self._annotate_row(rna_entry, dna_entry)
                    dna_entry = next(dna_reader, None)
                else:
                    yield rna_entry

    def annotate(self, stream):
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=[
            'Region',
            'Position',
            'Reference',
            'Strand',
            'Coverage',
            'MeanQ',
            'BaseCount[A,C,G,T]',
            'AllSubs',
            'Frequency',
            'gCoverage',
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
    parser.add_argument(
        '-d',
        '--debug',
        help='Report stack trace on crash.',
        action='store_true',
    )
    return parser.parse_args()


def main():
    options = parse_options()
    try:
        x = RTAnnotater(options.rna_file, options.dna_file)
        x.annotate(sys.stdout)
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
        sys.exit(1)
