import sys
import unittest
from contextlib import contextmanager
from io import StringIO
from reditools.tools.analyze.parse_args import parse_args

__all__ = ('TestParseArgs')


@contextmanager
def capture_sys_output():
    capture_out, capture_err = StringIO(), StringIO()
    current_out, current_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = capture_out, capture_err
        yield capture_out, capture_err
    finally:
        sys.stdout, sys.stderr = current_out, current_err


class TestParseArgs(unittest.TestCase):
    def test_legacy_pruning(self):
        args = parse_args(['test/test.bam'])
        self.assertFalse(hasattr(args, 'dna'))
        self.assertFalse(hasattr(args, 'exclude_multis'))
        self.assertFalse(hasattr(args, 'strict'))
        self.assertFalse(hasattr(args, 'load_omopolymeric_file'))

    def test_dna_mode(self):
        args = parse_args([
            'test/test.bam',
            '--dna',
        ])
        self.assertEqual(args.strand, 0)
        self.assertFalse(hasattr(args, 'dna'))

    def test_exclude_multis(self):
        args = parse_args([
            'test/test.bam',
            '--exclude-multis',
        ])
        self.assertEqual(args.max_editing_nucleotides, 1)
        self.assertFalse(hasattr(args, 'exclude_multis'))

    def test_strict(self):
        args = parse_args([
            'test/test.bam',
            '--strict',
        ])
        self.assertEqual(args.min_edits, 1)
        self.assertFalse(hasattr(args, 'strict'))

    def test_load_omopolymeric_file(self):
        args = parse_args([
            'test/test.bam',
            '--load-omopolymeric-file',
            'test/test.bed',
        ])
        self.assertEqual(args.exclude_regions, ['test/test.bed'])
        self.assertFalse(hasattr(args, 'load_omopolymeric_file'))

        args = parse_args([
            'test/test.bam',
            '--exclude-regions', 'example.bed',
            '--load-omopolymeric-file', 'test/test.bed',
        ])
        self.assertEqual(
            args.exclude_regions,
            ['example.bed', 'test/test.bed'],
        )
        self.assertFalse(hasattr(args, 'load_omopolymeric_file'))

    def test_edit_frequency(self):
        with self.assertRaises(SystemExit):
            with capture_sys_output() as (stdout, stderr):
                parse_args([
                    'test/test.bam',
                    '--max-editing-nucleotides', '1',
                    '--min-edits', '3',
                ])

    def test_unstranded(self):
        with self.assertRaises(SystemExit):
            with capture_sys_output() as (stdout, stderr):
                parse_args([
                    'test/test.bam',
                    '--strand', '0',
                    '--strand-correction',
                ])
