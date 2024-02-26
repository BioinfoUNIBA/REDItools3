"""Miscellaneous utility functions."""

import csv
import os
from gzip import open as gzip_open


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

    Returns:
        Iterator of BED file contents.
    """
    stream = open_stream(path)
    return csv.reader(stream, delimiter='\t')


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
