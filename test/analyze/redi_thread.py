"""Test cases for setup_rtools() method."""
from __future__ import annotations

import pathlib
import unittest
from test.sam_gen import SAM, Sequence, ntf

from reditools.logger import Logger
from reditools.reditools import REDItools
from reditools.tools.analyze.parse_args.parse_args import parse_args
from reditools.tools.analyze.redi_thread import (
    REDIThread,
    REDIThreadDNA,
    contig_order_from_bam,
    log_level,
)


class TestREDIThread(unittest.TestCase):
    """Test cases for setup_rtools() method."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.genome_fname = ntf(suffix=".fa")
        self.bam_fname = ntf(suffix=".bam")

        sam = SAM()
        sam.add_contig("chr2", length=60)
        sam.add_read("chr2", Sequence(
            sam.genome["chr2"],
            0,
        ))
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
        """Check that CLI options properly setup the thread pieces."""
        options = parse_args([
            self.bam_fname,
            "-r", "genome.fa",
            "--min-base-position", "10",
            "--max-base-position", "30",
            "--min-base-quality", "23",
            "--strand", "1",
            "--strand-confidence-threshold", "0.567",
            "--strand-correction",
            "--min-read-length", "33",
            "--min-read-quality", "99",
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

        sam_manager = redi_thread.sam_manager
        self.assertEqual(sam_manager.min_quality, 99)
        self.assertEqual(sam_manager.min_length, 33)


    def test_log_level(self) -> None:
        """Check CLI options for log level."""
        options = parse_args([self.bam_fname])
        self.assertEqual(log_level(options), Logger.silent_level)
        options = parse_args([self.bam_fname, "--verbose"])
        self.assertEqual(log_level(options), Logger.info_level)
        options = parse_args([self.bam_fname, "--debug"])
        self.assertEqual(log_level(options), Logger.debug_level)

    def test_contig_order_from_bam(self) -> None:
        """Test contig_order_from_bam() function."""
        contig_order = contig_order_from_bam(self.bam_fname)
        self.assertEqual(contig_order, {
            "chr2": 1,
            "chr1": 2,
        })

    def test_dna_options(self) -> None:
        """Check that CLI options properly setup the RNA-DNA thread pieces."""
        options = parse_args([
            self.bam_fname,
            "--dna-sample", self.bam_fname,
            "-r", "genome.fa",
            "--min-base-position", "10",
            "--max-base-position", "30",
            "--min-base-quality", "23",
            "--strand", "1",
            "--strand-confidence-threshold", "0.567",
            "--strand-correction",
            "--min-read-length", "33",
            "--min-read-quality", "99",

            "--dna-min-base-position", "20",
            "--dna-max-base-position", "40",
            "--dna-min-base-quality", "43",
            "--dna-min-read-length", "53",
            "--dna-min-read-quality", "999",
        ])

        redi_thread = REDIThreadDNA(options)

        rtools = redi_thread.rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, "genome.fa")
        self.assertEqual(rtools.min_base_position, 10)
        self.assertEqual(rtools.max_base_position, 30)
        self.assertEqual(rtools.min_base_quality, 23)
        self.assertEqual(rtools.strand, 1)
        self.assertEqual(rtools.strand_confidence_threshold, 0.567)
        self.assertTrue(rtools._use_strand_correction)

        sam_manager = redi_thread.sam_manager
        self.assertEqual(sam_manager.min_quality, 99)
        self.assertEqual(sam_manager.min_length, 33)

        rtools = redi_thread.dna_rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, "genome.fa")
        self.assertEqual(rtools.min_base_position, 20)
        self.assertEqual(rtools.max_base_position, 40)
        self.assertEqual(rtools.min_base_quality, 43)
        self.assertEqual(rtools.strand, 0)

        sam_manager = redi_thread.dna_sam_manager
        self.assertEqual(sam_manager.min_quality, 999)
        self.assertEqual(sam_manager.min_length, 53)
