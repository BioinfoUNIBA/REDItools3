"""Test cases for CompiledPosition class."""
from __future__ import annotations

import unittest

from reditools.compiled_position import (
    CompiledPosition,
    _calc_edit_ratio,
    _strand_from_string,
    rtresult_from_compiled_position,
)


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
        """Check variant summary for RTResult."""
        self.cp.add_base(10, "+", "A")
        self.cp.add_base(10, "+", "C")
        for _ in range(3):
            self.cp.add_base(10, "+", "T")
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.variants, ["T", "C"])

    def test_strand_from_string(self) -> None:
        """Check the _strand_from_string function."""
        self.assertEqual(_strand_from_string("+"), "+")
        self.assertEqual(_strand_from_string("-"), "-")
        self.assertEqual(_strand_from_string("*"), "*")
        self.assertEqual(_strand_from_string("1"), "+")
        self.assertEqual(_strand_from_string("0"), "-")
        self.assertEqual(_strand_from_string("2"), "*")

    def test_calc_edit_ratio(self) -> None:
        """Check the _calc_edit_ratio functions."""
        self.assertEqual(
            _calc_edit_ratio({"A": 0, "C": 0, "G": 0, "T": 10}, "T"),
            0,
        )
        self.assertEqual(
            _calc_edit_ratio({"A": 0, "C": 2, "G": 0, "T": 6}, "T"),
            0.25,
        )
        self.assertEqual(
            _calc_edit_ratio({"A": 0, "C": 2, "G": 1, "T": 6}, "T"),
            0.25,
        )

    def test_edit_ratio(self) -> None:
        """Check edit_ratio cacluation from RTResult."""
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0)

        self.cp.add_base(10, "+", "A")
        self.cp.add_base(10, "+", "C")
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "T")
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.5)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(10, "+", "C")
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.edit_ratio, 0.75)

    def test_mean_quality(self) -> None:
        """Check mean_quality calculation from RTResult."""
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 0)

        self.cp.add_base(10, "+", "C")
        self.cp.add_base(20, "+", "C")
        self.cp.add_base(30, "+", "C")
        rtresult = rtresult_from_compiled_position(self.cp, "*")
        self.assertEqual(rtresult.mean_quality, 20)

    def test_noise_threshold(self) -> None:
        """Check noise threshold feature."""
        for base in ("A", "A", "A", "T", "G", "G"):
            self.cp.add_base(10, "+", base)

        self.assertEqual(
            rtresult_from_compiled_position(self.cp, "*", 0.1).variants,
            ["G", "T"],
        )
        self.assertEqual(
            rtresult_from_compiled_position(self.cp, "*", 0.3).variants,
            ["G"],
        )
        self.assertEqual(
            rtresult_from_compiled_position(self.cp, "*", 0.4).variants,
            [],
        )

    def test_split_by_strand(self) -> None:
        """Check split_by_strand function."""
        self.cp.add_base(40, "+", "A")
        self.cp.add_base(35, "-", "C")
        self.cp.add_base(30, "+", "G")
        fwd, rev = self.cp.split_by_strand()
        self.assertEqual(fwd.bases, ["A", "G"])
        self.assertEqual(fwd.strands, ["+", "+"])
        self.assertEqual(rev.bases, ["C"])
        self.assertEqual(rev.strands, ["-"])

    def test_split_by_strand_no_reads(self) -> None:
        """Check split_by_strand function when there are no reads."""
        cps = self.cp.split_by_strand()
        self.assertEqual(len(cps), 0)

    def test_split_by_strand_one_stranded(self) -> None:
        """Check split_by_strand function when only one strand present."""
        self.cp.add_base(40, "-", "A")
        self.cp.add_base(40, "-", "A")
        self.cp.add_base(40, "-", "G")
        cps = self.cp.split_by_strand()
        self.assertEqual(len(cps), 1)
        self.assertEqual(cps[0].bases, ["A", "A", "G"])
        self.assertEqual(cps[0].strands, ["-", "-", "-"])
