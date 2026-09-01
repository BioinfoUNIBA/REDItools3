"""Test cases for setup_rtools() method."""
from __future__ import annotations

import pathlib
import unittest
import pysam
from test.sam_gen import SAM, Sequence, ntf

from reditools.logger import Logger
from reditools.reditools import REDItools
from reditools.region import Region
from reditools.tools.analyze.parse_args.parse_args import parse_args
from reditools.tools.analyze import redi_thread


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
        sam.add_contig("chr3", length=60)

        sam.genome.save_to_fasta(self.genome_fname)
        sam.save_to_sam(self.bam_fname, self.genome_fname)

        self.cli_input = [
            self.bam_fname,
            "-r", self.genome_fname,
            "--min-base-position", "10",
            "--max-base-position", "30",
            "--min-base-quality", "23",
            "--strand", "1",
            "--strand-confidence-threshold", "0.567",
            "--strand-correction",
            "--min-read-length", "33",
            "--min-read-quality", "99",
        ]

    def tearDown(self) -> None:
        """Posttest teardown."""
        pathlib.Path(self.genome_fname).unlink()
        pathlib.Path(f"{self.genome_fname}.fai").unlink()
        pathlib.Path(self.bam_fname).unlink()
        pathlib.Path(f"{self.bam_fname}.bai").unlink()

    def test_options(self) -> None:
        """Check that CLI options properly setup the thread pieces."""
        options = parse_args(self.cli_input)
        thread = redi_thread.REDIThread(options)
        rtools = thread.rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, self.genome_fname)
        self.assertEqual(rtools.min_base_position, 10)
        self.assertEqual(rtools.max_base_position, 30)
        self.assertEqual(rtools.min_base_quality, 23)
        self.assertEqual(rtools.strand, 1)
        self.assertEqual(rtools.strand_confidence_threshold, 0.567)
        self.assertTrue(rtools._use_strand_correction)

        sam_manager = thread.sam_manager
        self.assertEqual(sam_manager.min_quality, 99)
        self.assertEqual(sam_manager.min_length, 33)

    def test_log_level(self) -> None:
        """Check CLI options for log level."""
        options = parse_args([self.bam_fname])
        self.assertEqual(redi_thread.log_level(options), Logger.silent_level)
        options = parse_args([self.bam_fname, "--verbose"])
        self.assertEqual(redi_thread.log_level(options), Logger.info_level)
        options = parse_args([self.bam_fname, "--debug"])
        self.assertEqual(redi_thread.log_level(options), Logger.debug_level)

    def test_contig_order_from_bam(self) -> None:
        """Test contig_order_from_bam() function."""
        contig_order = redi_thread.contig_order_from_bam(self.bam_fname)
        self.assertEqual(contig_order, {
            "chr2": 1,
            "chr1": 2,
            "chr3": 3,
        })

    def test_dna_options(self) -> None:
        """Check that CLI options properly setup the RNA-DNA thread pieces."""
        options = parse_args([
            *self.cli_input,
            "--dna-sample", self.bam_fname,
            "--dna-min-base-position", "20",
            "--dna-max-base-position", "40",
            "--dna-min-base-quality", "43",
            "--dna-min-read-length", "53",
            "--dna-min-read-quality", "999",
        ])

        thread = redi_thread.REDIThreadDNA(options)

        rtools = thread.rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, self.genome_fname)
        self.assertEqual(rtools.min_base_position, 10)
        self.assertEqual(rtools.max_base_position, 30)
        self.assertEqual(rtools.min_base_quality, 23)
        self.assertEqual(rtools.strand, 1)
        self.assertEqual(rtools.strand_confidence_threshold, 0.567)
        self.assertTrue(rtools._use_strand_correction)

        sam_manager = thread.sam_manager
        self.assertEqual(sam_manager.min_quality, 99)
        self.assertEqual(sam_manager.min_length, 33)

        rtools = thread.dna_rtools
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, self.genome_fname)
        self.assertEqual(rtools.min_base_position, 20)
        self.assertEqual(rtools.max_base_position, 40)
        self.assertEqual(rtools.min_base_quality, 43)
        self.assertEqual(rtools.strand, 0)

        sam_manager = thread.dna_sam_manager
        self.assertEqual(sam_manager.min_quality, 999)
        self.assertEqual(sam_manager.min_length, 53)

    def test_thread_manager_init(self) -> None:
        """Thread manager should crash if not initalized."""
        # Have to do this hard-coded because the thread can be
        # initialized from other test code.
        redi_thread.REDIThreadManager.thread = None
        with self.assertRaises(redi_thread.UninitializedError):
            redi_thread.REDIThreadManager.analyze(
                Region("chr1", 0, 100),
                "example.out",
            )

    def test_thread_manager_done(self) -> None:
        """Check that the thread manager quites early if previously done."""
        result_file = ntf()
        done_file = f"{result_file}.done"
        pathlib.Path(done_file).touch()
        redi_thread.REDIThreadManager.init_thread(parse_args(self.cli_input))
        self.assertFalse(redi_thread.REDIThreadManager.analyze(
            Region("chr1", 0, 100),
            result_file,
        ))
        pathlib.Path(result_file).unlink()
        pathlib.Path(done_file).unlink()

    def test_dna_filter(self) -> None:
        """Check the filter option will only include annotated rows."""
        dna_sam_file = ntf(suffix=".bam")
        with pysam.FastaFile(self.genome_fname) as fa:
            sam = SAM()
            for _ in fa.references:
                sam.add_contig(_, sequence=fa.fetch(_))
            sam.add_read(
                "chr1",
                Sequence(sam.genome["chr1"][:5], 0),
            )
            sam.save_to_sam(dna_sam_file, self.genome_fname)

        thread = redi_thread.REDIThreadDNA(parse_args([
            self.bam_fname,
            "-r", self.genome_fname,
            "--dna-sample", dna_sam_file,
            "--dna-filter",
            "--min-edits", "0",
            "--dna-min-read-length", "0",
        ]))

        results_file = ntf()
        thread.analyze(Region("chr1", 0, 100), results_file)
        with pathlib.Path(results_file).open("r") as stream:
            self.assertEqual(len(list(stream)), 5)
        pathlib.Path(results_file).unlink()

        results_file = ntf()
        thread.analyze(Region("chr3", 0, 100), results_file)
        with pathlib.Path(results_file).open("r") as stream:
            self.assertEqual(len(list(stream)), 0)
        pathlib.Path(results_file).unlink()

        pathlib.Path(dna_sam_file).unlink()
        pathlib.Path(f"{dna_sam_file}.bai").unlink()
