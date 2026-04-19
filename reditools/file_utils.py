"""Miscellaneous utility functions."""

import csv
import os
from gzip import open as gzip_open
from typing import Iterator, IO
from reditools.region import Region

def open_stream(path: str, mode: str='rt', encoding: str='utf-8'):
    """
    Open a input or output stream from a file, accounting for gzip.

    Parameters:
        path (str): Path to file for reading or writing
        mode (str): File mode
        encoding (str): File encoding

    Returns:
        TextIOWrapper to the file
    """
    if path.endswith('gz'):
        return gzip_open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)  # noqa: WPS515


def read_bed_file(*path: str) -> Iterator[Region]:
    """
    Return an iterator for a BED file.

    Parameters:
        path (str): Path to a BED file for reading.

    Yields:
        BED file contents as Regions.
    """
    if len(path) > 1:
        yield from read_bed_file(*path[1:])
    with open_stream(path[0]) as stream:
        reader = csv.reader(
            filter(lambda row: row[0] != '#', stream),
            delimiter='\t',
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
        encoding: str='utf-8',
) -> None:
    """
    Combine one or more files into another file.

    Parameters:
        output (file): A file like object for writing
        *fnames (string): Paths to files for concatenation
        clean_up (bool): If True, deletes the files after concatenation
        encoding (string): File encoding
    """
    for fname in fnames:
        with open(fname, 'r', encoding=encoding) as stream:
            for line in stream:
                output.write(line)
        if clean_up:
            os.remove(fname)


def load_text_file(file_name: str) -> list[str]:
    """
    Extract file contents to a list.

    Parameters:
        file_name (str): The file to open.

    Returns:
        List of content
    """
    with open_stream(file_name, 'r') as stream:
        return [line.strip() for line in stream]


def _read_splice_sites(stream: IO) -> Iterator[tuple[str, int, str, str]]:
    reader = csv.reader(stream, delimiter=' ')
    for idx, row in enumerate(reader, start=1):
        if row[0].startswith('#'):
            continue
        try:  # noqa: WPS229
            assert len(row) == 5
            assert row[3] in ('A', 'D')
            assert row[4] in ('+', '-')
            position = int(row[1])
            yield (row[0], position, row[3], row[4])
        except (AssertionError, ValueError) as exc:
            raise ValueError(
                f'Cannot parse splice file entry ({stream.name}:{idx})'
            ) from exc

def _splice_site_to_region(
        contig: str,
        position: int,
        splice: str,
        strand: str,
        splicing_span: int,
) -> Region | None:
    strand_map = {'-': 'D', '+': 'A'}
    position = position - 1
    if strand_map[strand] == splice:
        start = max(position - splicing_span, 0)
        stop = position
    else:
        start = position
        stop = position + splicing_span
    if start != stop:
        return Region(contig=contig, start=start, stop=stop)
    return None

def load_splicing_file(
        splicing_file: str,
        splicing_span: int,
) -> Iterator[Region]:
    """
    Read splicing positions from a file.

    Parameters:
        splicing_file (str): File path
        splicing_span(int): Width of splice sites

    Yeilds:
        Splicing file contents as Regions.
    """

    with open_stream(splicing_file) as stream:
        for splice_data in _read_splice_sites(stream):
            region = _splice_site_to_region(*splice_data, splicing_span)
            if region is not None:
                yield region
