import os
import unittest
from tempfile import NamedTemporaryFile

from reditools.rtannotater import RTAnnotater


class TestRTAnnotater(unittest.TestCase):
    def test_legacy_translate(self):
        test_dict = {
            'Coverage-q30': '100',
            'gCoverage-q30': '200',
            'AnotherField': '123',
        }
        RTAnnotater.legacy_translate(test_dict)
        self.assertEqual(test_dict, {
            'Coverage': '100',
            'gCoverage': '200',
            'AnotherField': '123',
        })
   
    def test_cmp_position(self):
        contig_order = {
            'chrZ': 1,
            'chr1': 2,
            'chr2': 3,     
        }

        rta = RTAnnotater(contig_order)
        self.assertEqual(
            rta.cmp_position(
                {'Region': 'chr1', 'Position': '5'},
                {'Region': 'chr1', 'Position': '5'},
            ),
            0,
        )
        self.assertTrue(
            rta.cmp_position(
                {'Region': 'chr1', 'Position': '5'},
                None,
            ) < 0,
        )
        self.assertTrue(
            rta.cmp_position(
                {'Region': 'chrZ', 'Position': '1'},
                {'Region': 'chr2', 'Position': '1'},
            ) < 0,
        )
        self.assertTrue(
            rta.cmp_position(
                {'Region': 'chr1', 'Position': '10'},
                {'Region': 'chr1', 'Position': '1'},
            ) > 0,
        )

    def test_annotate_row(self):
        rta = RTAnnotater({})
        self.assertEqual(
            rta.annotate_row(
                {
                    'gCoverage': '-',
                    'gMeanQ': '-',
                    'gBaseCount[A,C,G,T]': '-',
                    'gAllSubs': '-',
                    'gFrequency': '-',
                    'AnotherField': 'ABCD',
                },
                {
                    'Coverage': '100',
                    'MeanQ': '40',
                    'BaseCount[A,C,G,T]': '[1, 2, 3, 4]',
                    'AllSubs': 'AC AG AT',
                    'Frequency': '0.5',
                    'AnotherField': 'EFGH',
                    'YetAnotherField': 'IJKL',
                },
            ),
            {
                'gCoverage': '100',
                'gMeanQ': '40',
                'gBaseCount[A,C,G,T]': '[1, 2, 3, 4]',
                'gAllSubs': 'AC AG AT',
                'gFrequency': '0.5',
                'AnotherField': 'ABCD',
            },
        )
                
    def test_merge_files(self): 
        fieldnames = [
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
            'gFrequency',
        ]

        with NamedTemporaryFile(
                delete=False,
                suffix='.out',
                mode='w',
        ) as stream:
            rna_file = stream.name
            stream.write('\t'.join(fieldnames))
            stream.write('\n')
            stream.write('\t'.join([
                'chr1', '3', 'A', '*', '5', '30', '[5, 0, 0, 0]', '', '0.0',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')
            stream.write('\t'.join([
                'chr1', '5', 'A', '*', '5', '30', '[2, 0, 3, 0]', 'AG', '0.6',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')
            stream.write('\t'.join([
                'chr2', '7', 'G', '*', '5', '30', '[2, 0, 3, 0]', 'GA', '0.4',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')

        with NamedTemporaryFile(
                delete=False,
                suffix='.out',
                mode='w',
        ) as stream:
            dna_file = stream.name
            stream.write('\t'.join(fieldnames))
            stream.write('\n')
            stream.write('\t'.join([
                'chrZ', '3', 'A', '*', '5', '30', '[5, 0, 0, 0]', '', '0.0',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')
            stream.write('\t'.join([
                'chr1', '5', 'A', '*', '5', '35', '[2, 0, 0, 3]', 'AT', '0.6',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')
            stream.write('\t'.join([
                'chr1', '9', 'G', '*', '5', '30', '[2, 0, 3, 0]', 'GA', '0.4',
                '-', '-', '-', '-', '-',
            ]))
            stream.write('\n')

        rta = RTAnnotater({'chrZ': 1, 'chr1': 2, 'chr2': 3})
        annotated_data = list(rta.merge_files(rna_file, dna_file))
        self.assertEqual(
            annotated_data.pop(0),
            {
                'Region': 'chr1',
                'Position': '3',
                'Reference': 'A',
                'Strand': '*',
                'Coverage': '5',
                'MeanQ': '30',
                'BaseCount[A,C,G,T]': '[5, 0, 0, 0]',
                'AllSubs': '',
                'Frequency': '0.0',
                'gCoverage': '-',
                'gMeanQ': '-',
                'gBaseCount[A,C,G,T]': '-',
                'gAllSubs': '-',
                'gFrequency': '-',
            },
        )
        self.assertEqual(
            annotated_data.pop(0),
            {
                'Region': 'chr1',
                'Position': '5',
                'Reference': 'A',
                'Strand': '*',
                'Coverage': '5',
                'MeanQ': '30',
                'BaseCount[A,C,G,T]': '[2, 0, 3, 0]',
                'AllSubs': 'AG',
                'Frequency': '0.6',
                'gCoverage': '5',
                'gMeanQ': '35',
                'gBaseCount[A,C,G,T]': '[2, 0, 0, 3]',
                'gAllSubs': 'AT',
                'gFrequency': '0.6',
            },
        )
        self.assertEqual(
            annotated_data.pop(0),
            {
                'Region': 'chr2',
                'Position': '7',
                'Reference': 'G',
                'Strand': '*',
                'Coverage': '5',
                'MeanQ': '30',
                'BaseCount[A,C,G,T]': '[2, 0, 3, 0]',
                'AllSubs': 'GA',
                'Frequency': '0.4',
                'gCoverage': '-',
                'gMeanQ': '-',
                'gBaseCount[A,C,G,T]': '-',
                'gAllSubs': '-',
                'gFrequency': '-',
            },
        )
        self.assertEqual(len(annotated_data), 0)
        os.remove(rna_file)
        os.remove(dna_file)
