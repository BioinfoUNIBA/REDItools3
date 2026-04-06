"""Miscellaneous utility functions."""

import csv
import os
from gzip import open as gzip_open

from reditools.region import Region

__all__ = (
    'open_stream',
    'read_bed_file',
    'concat',
    'load_splicing_file',
    'load_text_file',
)


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
    return open(path, mode, encoding=encoding)


def read_bed_file(path):
    """
    Return an iterator for a BED file.

    Parameters:
        path (str): Path to a BED file for reading.

    Yields:
        BED file contents as Regions.
    """
    with open_stream(path) as stream:
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
    reader = csv.reader(stream, delimiter=' ')
    for idx, row in enumerate(reader, start=1):
        if len(row) == 0 or row[0].startswith('#'):
            continue
        if len(row) != 5:
            raise ValueError(
                'Cannot parse splice site. Row must have 5 values '
                f'({stream.name}:{idx})'
            )
        try:
            position = int(row[1])
        except ValueError as exc:
            raise ValueError(
                f'Splice site must be an integer ({stream.name}:{idx})'
            ) from exc

        if row[3] not in ('A', 'D'):
            raise ValueError(
                f'Splice type must be A or D ({stream.name}:{idx})'
            )
        if row[4] not in ('+', '-'):
            raise ValueError(
                f'Strand must be + or - ({stream.name}:{idx})'
            )

        yield (row[0], position, row[3], row[4])


def load_splicing_file(splicing_file, splicing_span):
    """
    Read splicing positions from a file.

    Parameters:
        splicing_file (str): File path
        splicing_span(int): Width of splice sites

    Yeilds:
        Splicing file contents as Regions.
    """
    strand_map = {'-': 'D', '+': 'A'}

    with open_stream(splicing_file) as stream:
        for contig, position, splice, strand in _read_splice_sites(stream):
            position = position - 1
            if strand_map[strand] == splice:
                start = max(position - splicing_span, 0)
                stop = position
            else:
                start = position
                stop = position + splicing_span
            if start != stop:
                yield Region(contig=contig, start=start, stop=stop)
