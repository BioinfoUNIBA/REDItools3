from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, ntf

from reditools.region import Region


class TestRegion(unittest.TestCase):

    def test_str(self) -> None:
        self.assertEqual(str(Region("chr1", 100, 200)), "chr1:101-200")
        self.assertEqual(str(Region("chr1", 0, None)), "chr1")
        self.assertEqual(str(Region("chr1", 50, None)), "chr1:51")

    def test_even_split(self) -> None:
        region = Region("chr1", 0, 1000)
        windows = region.split(250)
        self.assertEqual(len(windows), 4)
        self.assertEqual(windows[0], Region("chr1", 0, 250))
        self.assertEqual(windows[-1], Region("chr1", 750, 1000))

    def test_uneven_split(self) -> None:
        region = Region("chr1", 0, 950)
        windows = region.split(300)
        self.assertEqual(len(windows), 4)
        self.assertEqual(windows[0], Region("chr1", 0, 300))
        self.assertEqual(windows[1], Region("chr1", 300, 600))
        self.assertEqual(windows[2], Region("chr1", 600, 900))
        self.assertEqual(windows[3], Region("chr1", 900, 950))

    def test_impossible_split(self) -> None:
        region = Region("chr1", 0, 100)
        windows = region.split(200)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0], Region("chr1", 0, 100))

    def test_nonzero_split(self) -> None:
        region = Region("chr2", 5, 122)
        windows = region.split(50)
        self.assertEqual(len(windows), 3)
        self.assertEqual(windows[0], Region("chr2", 5, 55))
        self.assertEqual(windows[1], Region("chr2", 55, 105))
        self.assertEqual(windows[2], Region("chr2", 105, 122))

    def test_none_split(self) -> None:
        with self.assertRaises(IndexError):
            Region("chr1", None, 100).split(50)
        with self.assertRaises(IndexError):
            Region("chr1", 50, None).split(50)

    def test_from_string(self) -> None:
        fasta_fname = ntf(suffix=".fa")
        bam_fname = ntf(suffix=".bam")

        sam_obj = SAM()
        chr1_len = 600
        sam_obj.add_contig("chr1", length=chr1_len)

        sam_obj.genome.save_to_fasta(fasta_fname)
        sam_obj.save_to_sam(bam_fname, fasta_fname)

        region = Region.from_string("chr1:101-200", bam_fname)
        self.assertEqual(region, Region("chr1", 100, 200))
        region = Region.from_string("chr1:104", bam_fname)
        self.assertEqual(region, Region("chr1", 103, chr1_len))
        region = Region.from_string("chr1", bam_fname)
        self.assertEqual(region, Region("chr1", 0, chr1_len))

        Path(fasta_fname).unlink()
        Path(bam_fname).unlink()

    def test_parse_string(self) -> None:
        region = Region.parse_string("chr1:101-200")
        self.assertEqual(region, ("chr1", 100, 200))

        with self.assertRaises(ValueError):
            Region.parse_string("chr1:-2")

    def test_to_int(self) -> None:
        self.assertEqual(Region._to_int("10"), 10)
        self.assertEqual(Region._to_int("10,000"), 10000)
        with self.assertRaises(ValueError):
            Region._to_int("X")

    def test_order(self) -> None:
        regions_list = [
            Region("chr1", 20, 30),
            Region("chr1", 10, 30),
            Region("chr1", 10, 20),
        ]
        self.assertEqual(sorted(regions_list), list(reversed(regions_list)))
