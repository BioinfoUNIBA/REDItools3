"""Miscellaneous utility functions."""

import csv
import os
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
    return csv.reader(
        filter(lambda row: row[0] != '#', stream),
        delimiter='\t',
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

def load_poly_positions(fname):
    """
    Read omopolymeric positions from a file.

    Parameters:
        fname (str): File path
    """
    poly_positions = defaultDict(set)
    with read_bed_file(fname) as reader:
        for row in reader:
            poly_positions[row[0]] = Region(
                contig=row[0],
                start=row[1],
                stop=row[2],
            )
    return poly_positions


def load_splicing_file(self, splicing_file, span):
    """
    Read splicing positions from a file.

    Parameters:
        splicing_file (str): File path
        span(int): Width of splice sites
    """
    splice_positions = defaultdict(SortedSet)
    strand_map = {'-': 'D', '+': 'A'}

    with open_stream(splicing_file, 'r') as stream:
        for line in stream:
            fields = line.strip().split()

            chrom = fields[0]
            strand = fields[4]
            splice = fields[3]
            span = int(fields[1])

            total_array[chrom] += span

            coe = -1 if strand_map.get(strand, None) == splice else 1
            new_positions = [1 + span + coe * fctr for fctr in range(span)]
            splice_positions[chrom] |= new_positions
        return splice_positions

