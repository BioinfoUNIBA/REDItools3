import os
import unittest
from reditools.region import Region
from reditools.tools.analyze.region_args import region_args
from ..sam_gen import SAM
from tempfile import NamedTemporaryFile


class TestRegionArgs(unittest.TestCase):
    def setUp(self):
        with NamedTemporaryFile(delete=False, suffix='.fa') as f:
            self.fasta_fname = f.name
        with NamedTemporaryFile(delete=False, suffix='.bam') as f:
            self.bam_fname = f.name

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
        regions = region_args(self.bam_fname, None, 0)
        self.assertEqual(len(regions), 3)

    def test_region_input(self):
        input_region = Region(contig='chr1', start=1, stop=100)
        regions = region_args(self.bam_fname, input_region, 0)
        self.assertEqual(regions, [input_region])

    def test_region_window(self):
        input_region = Region(contig='chr1', start=1, stop=100)
        regions = region_args(self.bam_fname, input_region, 10)
        self.assertEqual(len(regions), 10)

    def test_bam_window(self):
        regions = region_args(self.bam_fname, None, 70)
        self.assertEqual(len(regions), 5)
