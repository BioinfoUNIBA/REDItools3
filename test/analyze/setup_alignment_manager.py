"""Test cases for setup_alignment_manager() method."""
from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, ntf

from reditools.tools.analyze.setup_alignment_manager import (
    setup_alignment_manager,
)


class TestSetupAlignmentManager(unittest.TestCase):
    """Test cases for setup_alignment_manager() method."""

    def test_setup(self) -> None:
        """Check setup_alignment_manager() method."""
        fasta_fname = ntf(suffix=".fa")
        bam_fname = ntf(suffix=".bam")

        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=120)
        sam_obj.add_contig("chr2", length=80)
        sam_obj.add_contig("chr3", length=60)

        sam_obj.genome.save_to_fasta(fasta_fname)
        sam_obj.save_to_sam(bam_fname, fasta_fname)


        exclusions_fname = ntf(suffix=".bed")
        with Path(exclusions_fname).open("w") as stream:
            stream.write("bad_read")

        rtam = setup_alignment_manager(
            [bam_fname],
            min_read_quality=50,
            min_read_length=123,
            exclusions_file=exclusions_fname,
        )

        self.assertEqual(rtam.min_quality, 50)
        self.assertEqual(rtam.min_length, 123)
        self.assertIn("bad_read", rtam.excluded_read_names)

        Path(fasta_fname).unlink()
        Path(bam_fname).unlink()
        Path(exclusions_fname).unlink()
