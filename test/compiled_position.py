"""Test cases for CompiledPosition class."""
from __future__ import annotations

import unittest

from reditools.compiled_position import CompiledPosition


class TestCompiledPosition(unittest.TestCase):
    """Test cases for CompiledPosition class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.cp = CompiledPosition(ref="A", contig="chr1", position=100)

    def test_len(self) -> None:
        """Check len() function."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        self.assertEqual(len(self.cp), 3)

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

    def test_variants(self) -> None:
        """Check variants() function."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        self.assertEqual(
            self.cp.variants(),
            ["AC", "AG"],
        )
