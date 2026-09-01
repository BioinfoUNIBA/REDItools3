"""Test cases for the RTResult class."""
from __future__ import annotations

import unittest

from reditools.compiled_position import CompiledPosition, RTResult


class TestRTResult(unittest.TestCase):
    """Test cases for the RTResult class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.cp = CompiledPosition(ref="A", contig="chr1", position=100)

    def test_len(self) -> None:
        """Check depth property."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.depth, 3)

    def test_base_counts(self) -> None:
        """Check base_counts property."""
        self.cp.add_base(5, "+", "A")
        self.cp.add_base(5, "+", "A")
        self.cp.add_base(6, "-", "C")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.base_counts, {
            "A": 2,
            "C": 1,
            "G": 0,
            "T": 0,
        })

    def test_reference(self) -> None:
        """Check reference base propagates to RTResult class."""
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.reference, "A")

    def test_variants(self) -> None:
        """Check variant summary for RTResult."""
        self.cp.add_base(10, "+", "C")
        self.cp.add_base(10, "+", "A")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.variants, ["AC"])

    def test_edit_ratio(self) -> None:
        """Check edit_ratio cacluation from RTResult."""
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0)

        self.cp.add_base(10, "+", "A")
        self.cp.add_base(10, "+", "C")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "T")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(10, "+", "C")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.75)

    def test_mean_quality(self) -> None:
        """Check mean_quality calculation from RTResult."""
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 0)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(20, "+", "C")
        self.cp.add_base(30, "+", "C")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 20)

    def test_from_dict(self) -> None:
        """Check from_dict() method."""
        rtresult = RTResult.from_dict({
            "Region": "chr1",
            "Position": "137622",
            "Reference": "G",
            "Strand": "+",
            "Coverage": "1",
            "MeanQ": "38.00",
            "BaseCount[A,C,G,T]": "[1, 0, 0, 0]",
            "AllSubs": "GA",
            "Frequency": "1.00",
            "gCoverage": "-",
            "gMeanQ": "-",
            "gBaseCount[A,C,G,T]": "-",
            "gAllSubs": "-",
            "gFrequency": "-",
        })
        self.assertEqual(rtresult.contig, "chr1")
        self.assertEqual(rtresult.position, 137621)
        self.assertEqual(rtresult.reference, "G")
        self.assertEqual(rtresult.strand, "+")
        self.assertEqual(rtresult.depth, 1)
        self.assertEqual(rtresult.mean_quality, 38)
        self.assertEqual(
            rtresult.base_counts,
            {"A": 1, "C": 0, "G": 0, "T": 0},
        )
        self.assertEqual(rtresult.variants, ["GA"])
        self.assertEqual(rtresult.edit_ratio, 1)

    def test_from_dict_number_strands(self) -> None:
        """Test from_dict() when strands are numbers."""
        rtresult = RTResult.from_dict({
            "Region": "chr1",
            "Position": "137622",
            "Reference": "G",
            "Strand": "+",
            "Coverage": "1",
            "MeanQ": "38.00",
            "BaseCount[A,C,G,T]": "[1, 0, 0, 0]",
            "AllSubs": "GA",
            "Frequency": "1.00",
            "gCoverage": "-",
            "gMeanQ": "-",
            "gBaseCount[A,C,G,T]": "-",
            "gAllSubs": "-",
            "gFrequency": "-",
        })
        self.assertEqual(rtresult.strand, "+")

    def test_complement(self) -> None:
        """Test complement() method."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        rtresult.complement()
        self.assertEqual(rtresult.reference, "T")
        self.assertEqual(
            rtresult.base_counts,
            {"A": 0, "C": 1, "G": 1, "T": 1},
        )
        self.assertEqual(rtresult.variants, ["TC", "TG"])

    def test_output_list(self) -> None:
        """Test output_list() method."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.output_list(), [
            "chr1", "101", "A", "*", "3", "35.00", "[1, 1, 1, 0]", "AC AG",
            "0.50", "-", "-", "-", "-", "-",
        ])

    def test_output_list_number_strand(self) -> None:
        """Test output_list() method when using number strands."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        rtresult = RTResult.from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.output_list(True), [  # noqa: FBT003
            "chr1", "101", "A", "2", "3", "35.00", "[1, 1, 1, 0]", "AC AG",
            "0.50", "-", "-", "-", "-", "-",
        ])
        rtresult = RTResult.from_compiled_position(self.cp, "-")
        self.assertEqual(rtresult.output_list(True)[3], "0")  # noqa: FBT003
        rtresult = RTResult.from_compiled_position(self.cp, "+")
        self.assertEqual(rtresult.output_list(True)[3], "1")  # noqa: FBT003

    def test_annotate(self) -> None:
        """Test adding DNA result as annotation."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        rna = RTResult.from_compiled_position(self.cp, "*")

        dna_cp = CompiledPosition(ref="A", contig="chr1", position=100)
        dna_cp.add_base(40, "+", "A")
        dna_cp.add_base(38, "+", "A")
        dna = RTResult.from_compiled_position(dna_cp, "*")

        rna.dna = dna

        self.assertEqual(rna.output_list(), [
            "chr1", "101", "A", "*", "3", "35.00", "[1, 1, 1, 0]", "AC AG",
            "0.50", "2", "39.00", "[2, 0, 0, 0]", "-", "0.00",
        ])
