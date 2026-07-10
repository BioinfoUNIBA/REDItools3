"""File handling utilities."""
from __future__ import annotations

import csv
import tempfile
from gzip import open as gzip_open
from pathlib import Path
from typing import IO, Iterator

from reditools.region import Region


def open_stream(  # type: ignore[no-untyped-def] # noqa: ANN201
        path: str,
        mode: str="rt",
        encoding: str="utf-8",
):
    """Open a file stream, handling both plain and gzipped files.

    Parameters
    ----------
    path : str
        The path to the file.
    mode : str, optional
        The mode in which the file is opened (default is 'rt').
    encoding : str, optional
        The encoding to use (default is 'utf-8').

    Returns
    -------
    file-like object
        The opened file stream.
    """
    if path.endswith("gz"):
        return gzip_open(path, mode, encoding=encoding)
    return Path(path).open(mode, encoding=encoding)  # noqa: WPS515 SIM115

def read_bed_file(*path: str) -> Iterator[Region]:
    """Read genomic regions from one or more BED files.

    Parameters
    ----------
    *path : str
        Paths to the BED files.

    Yields
    ------
    Region
        The regions defined in the BED files.
    """
    if len(path) > 1:
        yield from read_bed_file(*path[1:])
    with open_stream(path[0]) as stream:
        reader = csv.reader(
            filter(lambda row: row[0] != "#", stream),
            delimiter="\t",
        )
        for row in reader:
            yield Region(
                contig=row[0],
                start=int(row[1]),
                stop=int(row[2]),
            )


def concat(
        output: IO,
        *fnames: str,
        clean_up: bool=True,
        encoding: str="utf-8",
) -> None:
    """Concatenate multiple files into a single output stream.

    Parameters
    ----------
    output : IO
        The output stream to write to.
    *fnames : str
        The names of the files to concatenate.
    clean_up : bool, optional
        Whether to remove the source files after concatenation
        (default is True).
    encoding : str, optional
        The encoding to use when reading files (default is 'utf-8').
    """
    for fname in fnames:
        with Path(fname).open("r", encoding=encoding) as stream:
            output.writelines(stream)
        if clean_up:
            Path(fname).unlink()


def load_text_file(file_name: str) -> list[str]:
    """Load lines from a text file into a list, stripping whitespace.

    Parameters
    ----------
    file_name : str
        The name of the file to load.

    Returns
    -------
    list[str]
        A list of stripped lines from the file.
    """
    with open_stream(file_name, "r") as stream:
        return [line.strip() for line in stream]

def make_dir(prefix: str | None=None, dirname: str | None=None) -> str:
    """Create a folder.

    Parameters
    ----------
    prefix : str
        Filename prefix.
    dirname : str
        Path to folder parent.

    Returns
    -------
    str
        Path to the folder.
    """
    with tempfile.NamedTemporaryFile(prefix=prefix, dir=dirname) as stream:
        valid_name = stream.name
    Path(valid_name).mkdir()
    return valid_name

