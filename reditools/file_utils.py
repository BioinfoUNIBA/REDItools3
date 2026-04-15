"""Miscellaneous utility functions."""

import csv
import os
import re
from gzip import open as gzip_open

from reditools.region import Region


def open_stream(path, mode='rt', encoding='utf-8'):
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


def read_bed_file(*path):
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


def concat(output, *fnames, clean_up=True, encoding='utf-8'):
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


def load_text_file(file_name):
    """
    Extract file contents to a list.

    Parameters:
        file_name (str): The file to open.

    Returns:
        List of content
    """
    with open_stream(file_name, 'r') as stream:
        return [line.strip() for line in stream]


def _read_splice_sites(stream):
    pa = re.compile(r'(\S+) (\d+) \d+ (A|D) ([+-])')
    for idx, row in enumerate(stream, start=1):
        if row.startswith('#'):
            continue
        try:
            yield pa.fullmatch(row.strip()).groups()
        except Exception as exc:
            raise ValueError(
                f'Cannot parse splice file entry ({stream.name}:{idx})'
            ) from exc

def _splice_site_to_region(contig, position, splice, strand, splicing_span):
    strand_map = {'-': 'D', '+': 'A'}
    position = int(position) - 1
    if strand_map[strand] == splice:
        start = max(position - splicing_span, 0)
        stop = position
    else:
        start = position
        stop = position + splicing_span
    if start != stop:
        return Region(contig=contig, start=start, stop=stop)

def load_splicing_file(splicing_file, splicing_span):
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
