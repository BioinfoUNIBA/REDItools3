import unittest
from reditools.region import Region
from pysam import AlignmentFile
from reditools.analyze.region_args import region_args

class TestRegionArgs(unittest.TestCase):
    def test_no_input(self):
        regions = region_args('test/test.bam', None, 0)
        self.assertEqual(len(regions), 194)

    def test_region_input(self):
        input_region = Region(contig='chr1', start=1, stop=100)
        regions = region_args('test/test.bam', input_region, 0)
        self.assertEqual(regions, [input_region])

    def test_region_window(self):
        input_region = Region(contig='chr1', start=1, stop=100)
        regions = region_args('test/test.bam', input_region, 10)
        self.assertEqual(len(regions), 10)

    def test_bam_window(self):
        regions = region_args('test/test.bam', None, 100000000)
        self.assertEqual(len(regions), 212)
