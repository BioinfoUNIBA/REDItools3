"""Test cases for splice_file module."""
from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

from reditools.region import Region
from reditools.splicing_file import load_splicing_file


class TestSplicingFile(unittest.TestCase):
    """Test cases for splice_file module."""

    def write_file(self, data_list: str | list, sep: str=" ") -> str:
        """Write data to a file.

        Parameters
        ----------
        data_list : str | list
            Data to write to file. If a list, will use sep to concatenate.
        sep : str
            Field seperator. Only applicable if data_list is a list.

        Returns
        -------
        str
            Path to output file.
        """
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8",
        ) as stream:
            for row in data_list:
                if isinstance(row, str):
                    stream.write(row)
                else:
                    stream.write(sep.join([str(_) for _ in row]))
                stream.write("\n")
            return stream.name

    def check_test_data(
        self,
        test_data: Iterable[list | tuple],
        real_data: list,
    ) -> None:
        """Perform consistency check between real and expected output.

        Parameters
        ----------
        test_data : Iterable[Iterable]
            Expected output. Uses the second sub element of each element.
        real_data : list
            Actual output.
        """
        self.assertEqual([_[1] for _ in test_data], real_data)

    def test_splicing_basic(self) -> None:
        """Check load_splicing_file() method."""
        test_data = [
            (
                ("chr1", "10", "25", "A", "+"),
                Region(contig="chr1", start=4, stop=9),
            ),
            (
                ("chr2", "20", "25", "D", "-"),
                Region(contig="chr2", start=14, stop=19),
            ),
            (
                ("chr3", "5", "15", "A", "-"),
                Region(contig="chr3", start=4, stop=9),
            ),
            (
                ("chr3", "5", "10", "D", "+"),
                Region(contig="chr3", start=4, stop=9),
            ),
        ]
        fname = self.write_file(
            ["#Header"] + [_[0] for _ in test_data],
        )
        splice_sites = list(load_splicing_file(fname, 5))
        self.check_test_data(test_data, splice_sites)
        Path(fname).unlink()

    def test_splicing_edge(self) -> None:
        """Check effects when at the very start of a contig."""
        test_data = [
            ("chr1", "1", "25", "A", "+"),
            ("chr1", "1", "25", "D", "-"),
            ("chr1", "3", "25", "D", "-"),
        ]
        fname = self.write_file(["#Header", *test_data])
        splice_sites = list(load_splicing_file(fname, 5))
        self.assertEqual(
            splice_sites,
            [Region(contig="chr1", start=0, stop=2)],
        )
        Path(fname).unlink()
