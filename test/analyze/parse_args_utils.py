"""Test cases for analyze bounded_types module."""
from __future__ import annotations

import unittest

from reditools.tools.analyze.parse_args.bounded_types import (
    CastFloatError,
    CastIntError,
    ValueAboveMaximumError,
    ValueBelowMinimumError,
    bounded_float,
    bounded_int,
    check_number_bounds,
)


class TestParseArgsUtils(unittest.TestCase):
    """Test cases for analyze bounded_types module."""

    def test_check_number_bounds_valid(self) -> None:
        """Check check_number_bounds() method for acceptable values."""
        check_number_bounds(5, min_value=1, max_value=10)
        check_number_bounds(1, min_value=1)
        check_number_bounds(10, max_value=10)

    def test_check_number_bounds_too_low(self) -> None:
        """Check check_number_bounds() with below minimum input."""
        with self.assertRaises(ValueBelowMinimumError):
            check_number_bounds(0, min_value=1)

    def test_check_number_bounds_too_high(self) -> None:
        """Check check_number_bounds() with above maximum input."""
        with self.assertRaises(ValueAboveMaximumError):
            check_number_bounds(11, max_value=10)

    def test_bounded_int_valid(self) -> None:
        """Check bounded_int() method with valid input."""
        conv = bounded_int(min_value=2, max_value=6)
        self.assertEqual(conv("4"), 4)
        self.assertEqual(conv("2"), 2)
        self.assertEqual(conv("6"), 6)

    def test_bounded_int_invalid_type(self) -> None:
        """Check bounded_int() method with non-int input."""
        conv = bounded_int()
        with self.assertRaises(CastIntError):
            conv("2.1")
        with self.assertRaises(CastIntError):
            conv("foo")

    def test_bounded_int_out_of_bounds(self) -> None:
        """Check bounded_int method with out of bounds input."""
        conv = bounded_int(min_value=3)
        with self.assertRaises(ValueBelowMinimumError):
            conv("1")
        conv = bounded_int(max_value=1)
        with self.assertRaises(ValueAboveMaximumError):
            conv("2")

    def test_bounded_float_valid(self) -> None:
        """Check bounded_float() method with valid input."""
        conv = bounded_float(min_value=0.5, max_value=2.6)
        self.assertEqual(conv("1.2"), 1.2)
        self.assertEqual(conv("0.5"), 0.5)
        self.assertEqual(conv("2.6"), 2.6)
        self.assertEqual(conv("1"), 1.0)

    def test_bounded_float_invalid_type(self) -> None:
        """Check bounded_float() with non-float input."""
        conv = bounded_float()
        with self.assertRaises(CastFloatError):
            conv("hello")

    def test_bounded_float_out_of_bounds(self) -> None:
        """Check bounded_float() with out of bounds input."""
        conv = bounded_float(min_value=0.1, max_value=1.2)
        with self.assertRaises(ValueBelowMinimumError):
            conv("0.01")
        with self.assertRaises(ValueAboveMaximumError):
            conv("2.0")
