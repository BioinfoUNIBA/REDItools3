import unittest
import argparse
from reditools.tools.analyze.parse_args import (
    check_number_bounds, bounded_int, bounded_float
)


class TestParseArgsUtils(unittest.TestCase):
    def test_check_number_bounds_valid(self):
        # No error for value in bounds
        check_number_bounds(5, min=1, max=10)
        check_number_bounds(1, min=1)
        check_number_bounds(10, max=10)

    def test_check_number_bounds_too_low(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            check_number_bounds(0, min=1)

    def test_check_number_bounds_too_high(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            check_number_bounds(11, max=10)

    def test_bounded_int_valid(self):
        conv = bounded_int(min=2, max=6)
        self.assertEqual(conv('4'), 4)
        self.assertEqual(conv('2'), 2)
        self.assertEqual(conv('6'), 6)

    def test_bounded_int_invalid_type(self):
        conv = bounded_int()
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('foo')

    def test_bounded_int_out_of_bounds(self):
        conv = bounded_int(min=3)
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('1')
        conv = bounded_int(max=1)
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('2')

    def test_bounded_float_valid(self):
        conv = bounded_float(min=0.5, max=2.6)
        self.assertEqual(conv('1.2'), 1.2)
        self.assertEqual(conv('0.5'), 0.5)
        self.assertEqual(conv('2.6'), 2.6)

    def test_bounded_float_invalid_type(self):
        conv = bounded_float()
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('hello')

    def test_bounded_float_out_of_bounds(self):
        conv = bounded_float(min=0.1, max=1.2)
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('0.01')
        with self.assertRaises(argparse.ArgumentTypeError):
            conv('2.0')
