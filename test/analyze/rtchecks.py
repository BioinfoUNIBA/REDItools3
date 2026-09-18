"""Test cases for RTChecks class."""
from __future__ import annotations

import unittest
from argparse import Namespace
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.compiled_position import CompiledPosition, RTResult
from reditools.tools.analyze.rtchecks import RTChecks


class TestRTChecks(unittest.TestCase):
    """Test cases for RTChecks class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.bases = CompiledPosition(contig="chr1", position=1, ref="A")
        self.options = Namespace(
            max_editing_nucleotides=4,
            min_read_depth=0,
            min_edits=0,
            min_edits_per_nucleotide=0,
            variants=["all"],
            exclude_regions=None,
            splicing_file=None,
            bed_file=None,
        )

    def run_check(self, is_none: bool) -> None:
        """Perform a quality control check.

        Create an RTChecks object from self.options and runs check(self.bases).
        Then checks whether the returned value is None.

        Parameters
        ----------
        is_none : bool
            If True, runs assertIsNone. Else, runs assesrtIsNotNone.
        """
        rtc = RTChecks(self.options)
        check_value =  rtc.check(RTResult(self.bases, "*"))
        if is_none:
            self.assertIsNone(check_value)
        else:
            self.assertIsNotNone(check_value)

    def test_check_column_edit_frequency(self) -> None:
        """Check --min-edits input."""
        self.options.min_edits = 1
        self.run_check(False)

        self.bases.add_base(quality=30, base="A", strand="*")
        self.run_check(False)

        self.options.min_edits = 0
        self.run_check(True)

        self.options.min_edits = 1
        for _ in range(4):
            self.bases.add_base(quality=30, base="T", strand="*")
        self.run_check(True)

        self.options.min_edits = 3
        self.run_check(False)

    def test_check_column_min_edits(self) -> None:
        """Check min-edits-per-nucleotide input."""
        self.options.min_edits_per_nucleotide = 1
        self.run_check(True)

        self.bases.add_base(quality=30, base="A", strand="*")
        self.bases.add_base(quality=30, base="A", strand="*")
        self.run_check(True)

        self.options.min_edits_per_nucleotide = 2
        self.run_check(True)

        self.bases.add_base(quality=30, base="T", strand="*")
        self.run_check(False)

        self.bases.add_base(quality=30, base="T", strand="*")
        self.run_check(True)

        self.bases.add_base(quality=30, base="C", strand="*")
        self.run_check(False)

    def test_check_min_read_depth(self) -> None:
        """Check min-read-depth input."""
        self.options.min_read_depth = 2
        self.run_check(False)

        self.bases.add_base(quality=30, base="C", strand="*")
        self.run_check(False)

        self.bases.add_base(quality=30, base="A", strand="*")
        self.run_check(True)

        self.bases.add_base(quality=30, base="A", strand="*")
        self.run_check(True)

    def test_check_exclusions(self) -> None:
        """Check exclude-regions input."""
        with NamedTemporaryFile(
                delete=False,
                suffix=".bed",
                mode="w+",
        ) as stream:
            stream.write("chr1\t20\t30\n")
            bed_file = stream.name

        self.options.exclude_regions = [bed_file]

        self.run_check(True)

        with Path(bed_file).open("a") as stream:
            stream.write("chr1\t0\t10\n")
        self.run_check(False)

        with Path(bed_file).open("a") as stream:
            stream.write("chr2\t0\t10\n")
        self.run_check(False)

        Path(bed_file).unlink()

    def test_check_splicing(self) -> None:
        """Check splicing-file input."""
        with NamedTemporaryFile(
            delete=False,
            suffix=".txt",
            mode="w+",
        ) as stream:
            stream.write("chr1 1 4 A +\n")
            splice_file = stream.name

        self.options.splicing_file = splice_file
        self.options.splicing_span = 4

        self.run_check(True)

        with Path(splice_file).open("w") as stream:
            stream.write("chr1 1 4 A -\n")
        self.run_check(False)

        with Path(splice_file).open("w") as stream:
            stream.write("chr1 1 4 D +\n")
        self.run_check(False)

        with Path(splice_file).open("w") as stream:
            stream.write("chr1 1 4 D -\n")
        self.run_check(True)

        Path(splice_file).unlink()

    def test_check_max_editing_nucleotides(self) -> None:
        """Check max-editing-nucleotides input."""
        self.options.max_editing_nucleotides = 1
        self.run_check(True)

        self.bases.add_base(quality=30, base="A", strand="*")
        self.run_check(True)

        self.bases.add_base(quality=30, base="T", strand="*")
        self.bases.add_base(quality=30, base="T", strand="*")
        self.run_check(True)

        self.bases.add_base(quality=30, base="C", strand="*")
        self.run_check(False)

        self.options.max_editing_nucleotides = 2
        self.run_check(True)

    def test_check_target_positions(self) -> None:
        """Check bed-file input."""
        with NamedTemporaryFile(
                delete=False,
                suffix=".bed",
                mode="w+",
        ) as stream:
            stream.write("chr1\t10\t20\n")
            bed_file = stream.name
        self.options.bed_file = [bed_file]
        self.run_check(False)

        with Path(bed_file).open("a") as stream:
            stream.write("chr1\t0\t20\n")
        self.run_check(True)

        with Path(bed_file).open("a") as stream:
            stream.write("chr2\t0\t20\n")
        self.run_check(True)

        Path(bed_file).unlink()
