"""Test cases for analyze parse_args."""
from __future__ import annotations

import sys
import unittest
from contextlib import contextmanager
from io import StringIO
from typing import Iterator

from reditools import reditools
from reditools.tools.analyze.parse_args.parse_args import parse_args


class TestParseArgs(unittest.TestCase):
    """Test cases for analyze parse_args."""

    def test_legacy_pruning(self) -> None:
        """Check that legacy options are removed from the Namespace."""
        args = parse_args(["test/test.bam"])
        self.assertFalse(hasattr(args, "dna"))
        self.assertFalse(hasattr(args, "exclude_multis"))
        self.assertFalse(hasattr(args, "strict"))
        self.assertFalse(hasattr(args, "load_omopolymeric_file"))

    def test_dna_mode(self) -> None:
        """Check that --dna sets --strand to zero."""
        args = parse_args([
            "test/test.bam",
            "--dna",
        ])
        self.assertEqual(args.strand, reditools.UNSTRANDED_MODE)
        self.assertFalse(hasattr(args, "dna"))

    def test_exclude_multis(self) -> None:
        """Check that --exclude-multis sets --max-editing-nucleotides to one."""
        args = parse_args([
            "test/test.bam",
            "--exclude-multis",
        ])
        self.assertEqual(args.max_editing_nucleotides, 1)
        self.assertFalse(hasattr(args, "exclude_multis"))

    def test_strict(self) -> None:
        """Check that --strict sets --min-edits to one."""
        args = parse_args([
            "test/test.bam",
            "--strict",
        ])
        self.assertEqual(args.min_edits, 1)
        self.assertFalse(hasattr(args, "strict"))

    def test_load_omopolymeric_file(self) -> None:
        """Check that --load-omopolymeric-file appends to exclude_regions."""
        args = parse_args([
            "test/test.bam",
            "--load-omopolymeric-file",
            "test/test.bed",
        ])
        self.assertEqual(args.exclude_regions, ["test/test.bed"])
        self.assertFalse(hasattr(args, "load_omopolymeric_file"))

        args = parse_args([
            "test/test.bam",
            "--exclude-regions", "example.bed",
            "--load-omopolymeric-file", "test/test.bed",
        ])
        self.assertEqual(
            args.exclude_regions,
            ["example.bed", "test/test.bed"],
        )
        self.assertFalse(hasattr(args, "load_omopolymeric_file"))

    def test_edit_frequency(self) -> None:
        """Check for error when min-edits is above max-editing-nucleotides."""
        with self.assertRaises(SystemExit), self.capture_sys_output():
            parse_args([
                "test/test.bam",
                "--max-editing-nucleotides", "1",
                "--min-edits", "3",
            ])

    def test_unstranded(self) -> None:
        """Check for error when strand is zero whiel using strand-correction."""
        with self.assertRaises(SystemExit), self.capture_sys_output():
            parse_args([
                "test/test.bam",
                "--strand", "0",
                "--strand-correction",
            ])

    @contextmanager
    def capture_sys_output(self) -> Iterator[tuple[StringIO, StringIO]]:
        """Capture standard and error output.

        This method is used to silence the output for parse_args().

        Returns
        -------
        Iterator[tuple[StringIO, StringIO]]
            Stdout and stderr, respectively.
        """
        capture_out, capture_err = StringIO(), StringIO()
        current_out, current_err = sys.stdout, sys.stderr
        try:  # noqa: WPS229
            sys.stdout = capture_out
            sys.stderr = capture_err
            yield capture_out, capture_err
        finally:
            sys.stdout = current_out
            sys.stderr = current_err

