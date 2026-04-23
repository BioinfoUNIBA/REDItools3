import csv
import os
import unittest
from tempfile import NamedTemporaryFile

from reditools.rtindexer import RTIndexer


class TestRTIndexer(unittest.TestCase):
    test_data = [
        {
            'Region': 'chr1',
            'Position': 1,
            'Reference': 'A',
            'BaseCount[A,C,G,T]': '[10, 0, 0, 0]',
        },
        {
            'Region': 'chr1',
            'Position': 2,
            'Reference': 'A',
            'BaseCount[A,C,G,T]': '[0, 0, 10, 0]',
        },
        {
            'Region': 'chr1',
            'Position': 3,
            'Reference': 'G',
            'BaseCount[A,C,G,T]': '[0, 10, 10, 0]',
        },
    ]

    def setUp(self):
        with NamedTemporaryFile(
                delete=False,
                suffix='.out',
                mode='wt',
        ) as stream:
            self.output_filename = stream.name
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    'Region',
                    'Position',
                    'Reference',
                    'BaseCount[A,C,G,T]',
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(self.test_data)

        with NamedTemporaryFile(
                delete=False,
                suffix='.bed',
                mode='wt',
        ) as stream:
            self.bed_filename = stream.name
            stream.write('chr1\t0\t2\n')

    def tearDown(self):
        os.remove(self.output_filename)
        os.remove(self.bed_filename)

    def test_baseline(self):
        rti = RTIndexer()
        rti.add_rt_output(self.output_filename)
        self.assertEqual(rti.calc_index(), {
            'A-C': 0,
            'A-T': 0,
            'A-G': 50,
            'C-A': 0,
            'C-T': 0,
            'C-G': 0,
            'G-A': 0,
            'G-C': 50,
            'G-T': 0,
            'T-A': 0,
            'T-C': 0,
            'T-G': 0,
        })

    def test_region(self):
        rti = RTIndexer(region=('chr1', 100, 200))
        self.assertFalse(rti.do_ignore({'Region': 'chr1', 'Position': '150'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr1', 'Position': '50'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr1', 'Position': '250'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr2', 'Position': '150'}))

        rti =RTIndexer(region=('chr1', 100, None))
        self.assertFalse(rti.do_ignore({'Region': 'chr1', 'Position': '150'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr1', 'Position': '50'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr2', 'Position': '150'}))

    def test_targets(self):
        rti = RTIndexer()
        rti.add_target_from_bed(self.bed_filename)
        self.assertFalse(rti.do_ignore({'Region': 'chr1', 'Position': '1'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr1', 'Position': '2'}))
        self.assertTrue(rti.do_ignore({'Region': 'chr2', 'Position': '1'}))

    def test_exclusions(self):
        rti = RTIndexer()
        rti.add_exclusions_from_bed(self.bed_filename)
        self.assertTrue(rti.do_ignore({'Region': 'chr1', 'Position': '1'}))
        self.assertFalse(rti.do_ignore({'Region': 'chr1', 'Position': '2'}))
        self.assertFalse(rti.do_ignore({'Region': 'chr2', 'Position': '1'}))
