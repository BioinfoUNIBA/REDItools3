import os
import unittest
from test.sam_gen import SAM, ntf

from reditools.region import Region
from reditools.tools.analyze.parse_args import parse_args
from reditools.tools.analyze.region_args import region_args


class TestRegionArgs(unittest.TestCase):
    def setUp(self):
        self.fasta_fname = ntf(suffix='.fa')
        self.bam_fname = ntf(suffix='.bam')

        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=120)
        sam_obj.add_contig('chr2', length=80)
        sam_obj.add_contig('chr3', length=60)

        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

    def tearDown(self):
        os.remove(self.fasta_fname)
        os.remove(self.bam_fname)

    def test_no_input(self):
        options = parse_args([self.bam_fname])
        regions = region_args(options)
        self.assertEqual(len(regions), 3)

    def test_region_input(self):
        options = parse_args([self.bam_fname, '--region', 'chr1:1-100'])
        regions = region_args(options)
        self.assertEqual(regions, [Region('chr1', 0, 100)])

    def test_region_window(self):
        options = parse_args([
            self.bam_fname,
            '--region',
            'chr1:1-100',
            '--window',
            '10',
        ])
        regions = region_args(options)
        self.assertEqual(len(regions), 10)

    def test_bam_window(self):
        options = parse_args([self.bam_fname, '--window', '70'])
        regions = region_args(options)
        self.assertEqual(len(regions), 5)
