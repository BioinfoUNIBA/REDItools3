import unittest
import os
from reditools.region import Region
from .sam_gen import SAM
from tempfile import NamedTemporaryFile


class TestRegion(unittest.TestCase):

    def test_str(self):
        self.assertEqual(str(Region('chr1', 100, 200)), 'chr1:101-200')
        self.assertEqual(str(Region('chr1', 0, None)), 'chr1')
        self.assertEqual(str(Region('chr1', 50, None)), 'chr1:51')

    def test_split(self):
        # Standard case: evenly divisible
        region = Region('chr1', 0, 1000)
        windows = region.split(250)
        self.assertEqual(len(windows), 4)
        self.assertEqual(windows[0], Region('chr1', 0, 250))
        self.assertEqual(windows[-1], Region('chr1', 750, 1000))

        # Edge case: region size not divisible by window size
        region = Region('chr1', 0, 950)
        windows = region.split(300)
        self.assertEqual(len(windows), 4)
        self.assertEqual(windows[0], Region('chr1', 0, 300))
        self.assertEqual(windows[1], Region('chr1', 300, 600))
        self.assertEqual(windows[2], Region('chr1', 600, 900))
        # last window smaller than others
        self.assertEqual(windows[3], Region('chr1', 900, 950))

        # Edge case: very small region, window size larger than region
        region = Region('chr1', 0, 100)
        windows = region.split(200)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0], Region('chr1', 0, 100))

        # Edge case: region starts at nonzero, not divisible
        region = Region('chr2', 5, 122)
        windows = region.split(50)
        self.assertEqual(len(windows), 3)
        self.assertEqual(windows[0], Region('chr2', 5, 55))
        self.assertEqual(windows[1], Region('chr2', 55, 105))
        self.assertEqual(windows[2], Region('chr2', 105, 122))

        with self.assertRaises(IndexError):
            Region('chr1', None, 100).split(50)
        with self.assertRaises(IndexError):
            Region('chr1', 50, None).split(50)

    def test_contains(self):
        region = Region('chr1', 100, 200)
        self.assertTrue(region.contains('chr1', 150))
        self.assertFalse(region.contains('chr1', 99))
        self.assertFalse(region.contains('chr1', 200))
        self.assertFalse(region.contains('chr2', 150))

    def test_from_string(self):
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.fa') as fasta:
            fasta_fname = fasta.name
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.bam') as bam:
            bam_fname = bam.name

        sam_obj = SAM()
        chr1_len = 600
        sam_obj.add_contig('chr1', length=chr1_len)

        sam_obj.genome.save_to_fasta(fasta_fname)
        sam_obj.save_to_sam(bam_fname, fasta_fname)

        region = Region.from_string('chr1:101-200', bam_fname)
        self.assertEqual(region, Region('chr1', 100, 200))
        region = Region.from_string('chr1:104', bam_fname)
        self.assertEqual(region, Region('chr1', 103, chr1_len))
        region = Region.from_string('chr1', bam_fname)
        self.assertEqual(region, Region('chr1', 0, chr1_len))

        os.remove(fasta_fname)
        os.remove(bam_fname)

    def test_parse_string(self):
        region = Region.parse_string('chr1:101-200')
        self.assertEqual(region, ('chr1', 100, 200))

        with self.assertRaises(ValueError):
            Region.parse_string('chr1:-2')

    def test_to_int(self):
        self.assertEqual(Region._to_int("10"), 10)
        self.assertEqual(Region._to_int("10,000"), 10000)
        with self.assertRaises(ValueError):
            Region._to_int("X")
