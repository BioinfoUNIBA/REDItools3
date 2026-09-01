"""Test cases for setup_rtools() method."""
from __future__ import annotations

import pathlib
import unittest
from test.sam_gen import SAM, Sequence, ntf

from reditools.reditools import REDItools
from reditools.tools.analyze.parse_args.parse_args import parse_args
from reditools.tools.analyze.redi_thread import REDIThread


class TestREDIThread(unittest.TestCase):
    """Test cases for setup_rtools() method."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.genome_fname = ntf(suffix=".fa")
        self.bam_fname = ntf(suffix=".bam")

        sam = SAM()
        sam.add_contig("chr1", length=60)
        sam.add_read("chr1", Sequence(
            sam.genome["chr1"],
            0,
        ))

        sam.genome.save_to_fasta(self.genome_fname)
        sam.save_to_sam(self.bam_fname, self.genome_fname)

    def tearDown(self) -> None:
        """Posttest teardown."""
        pathlib.Path(self.genome_fname).unlink()
        pathlib.Path(self.bam_fname).unlink()

    def test_options(self) -> None:
        """Check setup_rtools() method."""
        options = parse_args([
            self.bam_fname,
            "-r", "genome.fa",
            "--min-base-position", "10",
            "--max-base-position", "30",
            "--min-base-quality", "23",
            "--strand", "1",
            "--strand-confidence-threshold", "0.567",
            "--strand-correction",
        ])
        redi_thread = REDIThread(options)
        rtools = redi_thread.rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, "genome.fa")
        self.assertEqual(rtools.min_base_position, 10)
        self.assertEqual(rtools.max_base_position, 30)
        self.assertEqual(rtools.min_base_quality, 23)
        self.assertEqual(rtools.strand, 1)
        self.assertEqual(rtools.strand_confidence_threshold, 0.567)
        self.assertTrue(rtools._use_strand_correction)
