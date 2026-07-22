"""Test methods for analyze tool's write_results()."""
from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.compiled_position import CompiledPosition, RTResult
from reditools.logger import Logger
from reditools.tools.analyze.parse_args.parse_args import parse_args
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.write_results import write_results


class TestWriteResults(unittest.TestCase):
    """Test methods for analyze tool's write_results()."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.rtresults = []

        cp = CompiledPosition(ref="A", contig="chr1", position=100)
        cp.add_base(40, "+", "A")
        cp.add_base(30, "+", "G")
        cp.add_base(30, "+", "T")
        self.rtresults.append(RTResult.from_compiled_position(cp, "+"))
        self.rtresults.append(RTResult.from_compiled_position(cp, "-"))
        self.rtresults.append(RTResult.from_compiled_position(cp, "*"))

        with NamedTemporaryFile(mode="w", dir=".", delete=False) as temp_file:
            self.fname = temp_file.name
        self.log = Logger(Logger.silent_level).log

    def tearDown(self) -> None:
        """Post-checks cleanup."""
        Path(self.fname).unlink()

    def rtchecks(self, cli_args: list[str]) -> RTChecks:
        """Create RTChecks object from CLI arguments.

        Parameters
        ----------
        cli_args : list[str]
            Command line arguments.

        Returns
        -------
        RTChecks
            Checks for output filtering.
        """
        return RTChecks(parse_args(cli_args))

    def test_write_results(self) -> None:
        """Check basic functionality."""
        write_results(
            self.rtresults,
            self.fname,
            self.rtchecks(["x.bam"]),
            self.log,
        )
        with Path(self.fname).open("r") as stream:
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "+", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "-", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "*", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )

    def test_write_results_numbers(self) -> None:
        """Check strand-number option."""
        write_results(
            self.rtresults,
            self.fname,
            self.rtchecks(["x.bam"]),
            self.log,
            True,  # noqa: FBT003
        )
        with Path(self.fname).open("r") as stream:
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "1", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "0", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "2", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )

    def test_filter(self) -> None:
        """Check that RTChecks filters output."""
        write_results(
            self.rtresults,
            self.fname,
            self.rtchecks(["x.bam", "--min-read-depth", "5"]),
            self.log,
        )
        self.assertEqual(Path(self.fname).stat().st_size, 0)
