"""Test cases for RTAlignmentFile."""
from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, Sequence, ntf

from reditools.alignment_file import RTAlignmentFile


class TestRTAlignmentFile(unittest.TestCase):
    """Test cases for RTAlignmentFile."""

    def setUp(self) -> None:
        """Preflight setup."""
        self.sam_obj = SAM()
        self.sam_obj.add_contig("chr1", length=60)
        self.refseq = self.sam_obj.genome["chr1"]

        self.genome_fname = ntf(suffix=".fa")
        self.bam_fname = ntf(suffix=".bam")

    def tearDown(self) -> None:
        """Posttest teardown."""
        Path(self.genome_fname).unlink()
        Path(self.bam_fname).unlink()

    def test_fetch_by_position(self) -> None:
        """Check fetch_by_position standard functionality."""
        for start, stop in (
                (0, 20),
                (20, None),
                (20, None),
                (40, None),
        ):
            read_seq = self.refseq[start:stop]
            self.sam_obj.add_read("chr1", Sequence(read_seq, start))

        self.sam_obj.add_contig("chr2")
        self.sam_obj.add_read("chr2", Sequence(self.sam_obj.genome["chr2"], 0))

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads_iter = rtaf.fetch_by_position("chr1")
            self.assertEqual(len(next(reads_iter)), 1)
            self.assertEqual(len(next(reads_iter)), 2)
            self.assertEqual(len(next(reads_iter)), 1)

    def test_exclude_reads(self) -> None:
        """Check ability to exclude reads by name."""
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq, 0, read_name="exclude_me"),
        )
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq, 0, read_name="include_me"),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(
                self.bam_fname,
                excluded_read_names={"exclude_me"},
        ) as rtaf:
            reads = next(rtaf.fetch_by_position("chr1"))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].query_name, "include_me")

    def test_check_quality(self) -> None:
        """Check MAPQ quality filter."""
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq, 0, mapq=10),
        )
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq, 0, mapq=30, read_name="include_me"),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname, min_quality=20) as rtaf:
            reads = next(rtaf.fetch_by_position("chr1"))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].query_name, "include_me")

    def test_check_length(self) -> None:
        """Check minimum read length filter."""
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq[:20], 0),
        )
        self.sam_obj.add_read(
            "chr1",
            Sequence(self.refseq, 0, read_name="include_me"),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname, min_length=30) as rtaf:
            reads = next(rtaf.fetch_by_position("chr1"))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].query_name, "include_me")

    def test_check_se_flags(self) -> None:
        """Check for filtering by SAM flags."""
        for idx, flag in enumerate([0, 16]):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                read_name=f"se_good_{idx}",
            )
            self.sam_obj.add_read("chr1", read)
        for idx, flag in enumerate([4, 256, 272, 512, 1024, 2048, 2064]):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                read_name=f"se_bad_{idx}",
            )
            self.sam_obj.add_read("chr1", read)

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = next(rtaf.fetch_by_position("chr1"))
            self.assertEqual(len(reads), 2)
            for algn_seg in reads:
                self.assertTrue(
                    algn_seg.query_name is not None and \
                    algn_seg.query_name.startswith("se_good"),
                )

    def test_check_pe_flags(self) -> None:
        """Check for SAM paired end flags."""
        for idx, flag in enumerate([83, 99]):
            self.sam_obj.add_read_pair(
                "chr1",
                Sequence(
                    self.refseq,
                    0,
                    flag=flag,
                    read_name=f"pe_good_{idx}",
                ),
            )
        bad_flags = [73, 89, 137, 153, 329, 339, 345, 355, 393, 409]
        for idx, flag in enumerate(bad_flags):
            self.sam_obj.add_read_pair(
                "chr1",
                Sequence(
                    self.refseq,
                    0,
                    flag=flag,
                    read_name=f"pe_bad_{idx}",
                ),
            )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = next(rtaf.fetch_by_position("chr1"))
            self.assertEqual(len(reads), 4)
            for read in reads:
                self.assertTrue(
                    read.query_name is not None and \
                    read.query_name.startswith("pe_good"),
                )
