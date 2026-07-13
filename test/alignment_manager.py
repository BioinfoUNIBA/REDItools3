"""Test cases for AlignmentManager class."""
from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, Sequence, ntf

from reditools.alignment_manager import AlignmentManager


class TestAlignmentManager(unittest.TestCase):
    """Test cases for AlignmentManager class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.genome_fname = ntf(suffix=".fa")
        self.bam_fname1 = ntf(suffix=".bam")
        self.bam_fname2 = ntf(suffix=".bam")

        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=80)
        refseq = sam_obj.genome["chr1"]
        sam_obj.genome.save_to_fasta(self.genome_fname)

        sam_obj.add_read("chr1", Sequence(refseq, 0, read_name="1_1"))
        sam_obj.add_read("chr1", Sequence(refseq[20:], 20, read_name="1_2"))
        sam_obj.add_read("chr1", Sequence(refseq[40:], 40, read_name="1_3"))
        sam_obj.save_to_sam(self.bam_fname1, self.genome_fname)

        sam_obj = SAM()
        sam_obj.add_contig("chr1", sequence=refseq)
        sam_obj.add_read("chr1", Sequence(refseq[20:], 20, read_name="2_1"))
        sam_obj.add_read("chr1", Sequence(refseq[50:], 50, read_name="2_2"))
        sam_obj.save_to_sam(self.bam_fname2, self.genome_fname)

    def tearDown(self) -> None:
        """Post checks cleanup."""
        Path(self.genome_fname).unlink()
        Path(self.bam_fname1).unlink()
        Path(self.bam_fname2).unlink()

    def test_propagation(self) -> None:
        """Check that properties of AlignmentManager propagate to sub files."""
        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(self.bam_fname1)

        self.assertEqual(rtam._bams[0].readqc.min_length, 10)
        self.assertEqual(rtam._bams[0].readqc.min_quality, 30)

    def test_fetch_by_position(self) -> None:
        """Check fetch_by_position works for all sub files."""
        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(self.bam_fname1)
        rtam.add_file(self.bam_fname2)

        read_iter = rtam.fetch_by_position("chr1")

        read_group = next(read_iter)
        self.assertEqual(len(read_group), 1)
        self.assertEqual(read_group[0].query_name, "1_1")
        self.assertEqual(rtam.next_read_start, 20)

        read_group = next(read_iter)
        self.assertEqual(len(read_group), 2)
        self.assertIn("2_1", (_.query_name for _ in read_group))
        self.assertIn("1_2", (_.query_name for _ in read_group))
        self.assertEqual(rtam.next_read_start, 40)
