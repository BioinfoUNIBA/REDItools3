"""Test cases for REDItools class."""
from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, Sequence, ntf

from reditools import reditools
from reditools.alignment_manager import AlignmentManager
from reditools.compiled_position import CompiledPosition


class TestSplitREDItools(unittest.TestCase):
    """Test cases for REDItools class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.rtools = reditools.SplitREDItools()
        self.rtools.strand = reditools.FORWARD_STRAND_MODE

        self.bam_file = ntf(suffix=".bam")
        self.fa_file = ntf(suffix=".fa")

        self.sam_obj = SAM()
        self.sam_obj.add_contig("chr1", length=10)
        self.sam_obj.add_read("chr1", Sequence(self.sam_obj.genome["chr1"], 0))

        self.sam_obj.genome.save_to_fasta(self.fa_file)
        self.sam_obj.save_to_sam(self.bam_file, self.fa_file)

        self.rtam = AlignmentManager()
        self.rtam.add_file(self.bam_file)

    def tearDown(self) -> None:
        """Post-checks cleanup."""
        Path(self.bam_file).unlink()
        Path(self.fa_file).unlink()

    def test_process_bases(self) -> None:
        """Check _process_bases() function."""
        cp = CompiledPosition(ref="A", position=1, contig="chr1")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "T")
        cp.add_base(30, "+", "G")
        cp.add_base(30, "+", "G")

        fwd_result, rev_result = self.rtools._process_bases(cp)

        self.assertEqual(fwd_result.reference, "A")
        self.assertEqual(fwd_result.strand, "+")
        self.assertEqual(fwd_result.variants, ["AG"])

        self.assertEqual(rev_result.reference, "A")
        self.assertEqual(rev_result.strand, "-")
        self.assertEqual(rev_result.variants, ["AT"])

    def test_process_bases_monostrand(self) -> None:
        """Check _process_bases() function."""
        cp = CompiledPosition(ref="A", position=1, contig="chr1")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "A")
        cp.add_base(30, "-", "T")

        results_tuple = self.rtools._process_bases(cp)
        self.assertEqual(len(results_tuple), 1)
        rev_result = results_tuple[0]

        self.assertEqual(rev_result.reference, "A")
        self.assertEqual(rev_result.strand, "-")
        self.assertEqual(rev_result.variants, ["AT"])
