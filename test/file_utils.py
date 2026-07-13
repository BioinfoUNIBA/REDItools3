"""Test cases for file_utils."""
from __future__ import annotations

import gzip
import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

from reditools import file_utils
from reditools.region import Region


class TestFileUtils(unittest.TestCase):
    """Test cases for file_utils."""

    def write_file(self, data_list: str | Iterable, sep: str=" ") -> str:
        """Write data to file.

        Parameters
        ----------
        data_list : str | Iterable
            Data to write. If data_list is not a string, it will be joined
            using sep.
        sep : str
            Field seperator. Only used if data_list is not a string.

        Returns
        -------
        str
            Filename data saved to.
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
        test_data: Iterable[tuple],
        real_data: list,
    ) -> None:
        """Check test data matches real data.

        Parameters
        ----------
        test_data : Iterable[tuple]
            Assumes second element of each tuple what the output should be.
        real_data : list
            What the output actually was.
        """
        self.assertEqual([_[1] for _ in test_data], real_data)

    def test_open_stream_plain(self) -> None:
        """Check read/write plain text files."""
        test_str = "test123"
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8",
        ) as stream:
            stream.write(test_str)
            fname = stream.name
        with file_utils.open_stream(fname, "rt") as stream:
            file_content = stream.read()
        self.assertEqual(file_content, test_str)
        Path(fname).unlink()

    def test_open_stream_gzip(self) -> None:
        """Check read/write gzipped files."""
        test_str = "test_gzip"
        with NamedTemporaryFile(
                delete=False,
                suffix=".gz",
                mode="wb",
        ) as stream:
            stream.write(gzip.compress(bytes(test_str, "utf-8")))
            fname = stream.name
        with file_utils.open_stream(fname, "rt") as stream:
            file_content = stream.read()
        self.assertEqual(file_content, test_str)
        Path(fname).unlink()

    def test_read_bed_file(self) -> None:
        """Check read BED files."""
        bed_data = (
            (
                ("chr1", 10, 20),
                Region("chr1", 10, 20),
            ),
            (
                ("chr1", 30, 40),
                Region("chr1", 30, 40),
            ),
        )
        fname = self.write_file((_[0] for _ in bed_data), sep="\t")
        region_list = list(file_utils.read_bed_file(fname))
        self.check_test_data(bed_data, region_list)
        Path(fname).unlink()

    def test_read_many_bed_files(self) -> None:
        """Check read multiple BED files."""
        bed_data = (
            (
                ("chr1", 10, 20),
                Region("chr1", 10, 20),
            ),
            (
                ("chr1", 30, 40),
                Region("chr1", 30, 40),
            ),
        )
        fnames = [self.write_file([row[0]], sep="\t") for row in bed_data]
        region_list = list(file_utils.read_bed_file(*fnames))
        self.check_test_data(bed_data, sorted(region_list))
        for fname in fnames:
            Path(fname).unlink()

    def test_concat(self) -> None:
        """Check file concatenation."""
        file_contents = ("file1", "file2", "file3")
        file_names = [self.write_file([_]) for _ in file_contents]

        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8") as stream:
            file_utils.concat(stream, *file_names, encoding="utf-8")
            concat_filename = stream.name
        for fname in file_names:
            self.assertFalse(Path(fname).exists())

        with Path(concat_filename).open("r") as stream:
            self.assertEqual(
                stream.read(),
                "".join([f"{_}\n" for _ in file_contents]),
            )

        Path(concat_filename).unlink()

    def test_load_text_file(self) -> None:
        """Check read plaintext files."""
        text_lines = ["rowA", "rowB", "rowC"]
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8",
        ) as stream:
            fname = stream.name
            stream.write("\n".join(text_lines))
        loaded_text = file_utils.load_text_file(fname)
        self.assertEqual(loaded_text, text_lines)
        Path(fname).unlink()
