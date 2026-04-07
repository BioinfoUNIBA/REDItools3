from reditools import file_utils
import csv


class RTAnnotater:
    def __init__(self, rna_file, dna_file, contig_order):
        self.rna_file = rna_file
        self.dna_file = dna_file
        self.contig_order = contig_order

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
