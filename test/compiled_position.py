"""Test cases for CompiledPosition and RTResult classes."""
from __future__ import annotations

import unittest

from reditools.compiled_position import CompiledPosition, RTResult


class TestCompiledPosition(unittest.TestCase):
    """Test cases for CompiledPosition and RTResult classes."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.cp = CompiledPosition(ref="A", contig="chr1", position=100)

    def test_len(self) -> None:
        """Check len() functions for both classes."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        self.assertEqual(len(self.cp), 3)
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(len(rtresult), 3)

    def test_complement(self) -> None:
        """Check compelement() function of CompiledPosition."""
        self.cp.add_base(11, "+", "A")
        self.cp.add_base(12, "-", "C")
        self.cp.complement()
        self.assertEqual(self.cp.bases, ["T", "G"])
        self.assertEqual(self.cp.ref, "T")

    def test_calculate_strand(self) -> None:
        """Check calculate_strand() function of CompiledPosition."""
        self.cp.add_base(20, "+", "A")
        self.cp.add_base(21, "-", "A")
        self.cp.add_base(22, "+", "C")
        self.assertEqual(self.cp.calculate_strand(), "+")
        self.assertEqual(self.cp.calculate_strand(0.7), "*")

    def test_filter_by_strand(self) -> None:
        """Check fitler_by_strand() function."""
        self.cp.add_base(5, "+", "A")
        self.cp.add_base(5, "+", "A")
        self.cp.add_base(6, "-", "C")
        self.cp.filter_by_strand("+")
        self.assertEqual(len(self.cp), 2)
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult["A"], 2)
        self.assertEqual(rtresult["C"], 0)

    def test_filter_by_strand_star(self) -> None:
        """Check filter_by_strand() with undetermined strand."""
        self.cp.add_base(5, "*", "A")
        self.cp.add_base(5, "*", "A")
        self.cp.add_base(6, "+", "C")
        self.cp.add_base(6, "-", "T")
        self.assertEqual(
            self.cp.calculate_strand(threshold=1),
            "*",
        )
        self.cp.filter_by_strand("*")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(len(rtresult), 4)
        self.assertEqual(rtresult["A"], 2)
        self.assertEqual(rtresult["C"], 1)

    def test_reference(self) -> None:
        """Check reference base propagates to RTResult class."""
        self.assertEqual(self.cp.ref, "A")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.reference, "A")

    def test_get_base_counts(self) -> None:
        """Check RTResult base count summary."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "A")
        self.cp.add_base(30, "+", "C")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult["A"], 2)
        self.assertEqual(rtresult["C"], 1)
        self.assertEqual(rtresult["REF"], 2)

    def test_iter(self) -> None:
        """Check list casting for RTResult."""
        self.cp.add_base(41, "+", "A")
        self.cp.add_base(42, "+", "C")
        self.cp.add_base(43, "+", "G")
        counts = list(RTResult(self.cp, "*"))
        self.assertEqual(counts, [1, 1, 1, 0])

    def test_variants(self) -> None:
        """Check variant summary for RTResult."""
        self.cp.add_base(10, "+", "C")
        self.cp.add_base(10, "+", "A")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.variants, ["AC"])

    def test_edit_ratio(self) -> None:
        """Check edit_ratio cacluation from RTResult."""
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0)

        self.cp.add_base(10, "+", "A")
        self.cp.add_base(10, "+", "C")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "T")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(10, "+", "C")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.75)

    def test_mean_quality(self) -> None:
        """Check mean_quality calculation from RTResult."""
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 0)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(20, "+", "C")
        self.cp.add_base(30, "+", "C")
        rtresult = RTResult(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 20)
