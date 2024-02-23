"""Miscellaneous utility functions."""

import csv
import os
import re
import socket
from collections import defaultdict
from gzip import open as gzip_open

from pysam.libcalignmentfile import AlignmentFile
from sortedcontainers import SortedSet


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


def enumerate_positions(regions):
    """
    Convert a list of regions into a list of individual positions.

    Parameters:
        regions (list): A list of iterables. Each element must start
                        with a contig and start position. End position
                        is optional. Additional values will be ignored.

    Returns:
        Dictionary of contigs to SortedSets enumerating the individual
        positions.
    """
    positions = defaultdict(SortedSet)
    for reg in regions:
        contig = reg[0]
        start = int(reg[1]) - 1
        if len(reg) < 3:
            end = start + 1
        else:
            end = int(reg[2])
        positions[contig] |= range(start, end)
    return positions


def get_hostname_string():
    """
    Retrieve the machine hostname, ip, and proccess ID.

    Returns:
        String in the format "hostname|ip|pid"
    """
    hostname = socket.gethostname()
    ip_addr = socket.gethostbyname(hostname)
    pid = os.getpid()
    return f'{hostname}|{ip_addr}|{pid}'


def check_list(functions, **kwargs):
    """
    Run through a list of functions, determining if any return False.

    Parameters:
        functions (list): A list of function references
        **kwargs: Any arguments to be passed to the members of functions

    Returns:
        False if any function in check_list returns False, else True
    """
    for check in functions:
        if not check(**kwargs):
            return False
    return True


def to_int(string):
    """
    Convert a (potentially formatted) string to an int.

    Parameters:
        string (str): A string representation of an integer

    Returns:
        The integer values of the string.
    """
    return int(re.sub(r'[\s,]', '', string))


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


def get_contigs(sam_path):
    """
    Retrieve contig or chromsome data from an alignment file.

    Parameters:
        sam_path (string): Path to an alignment file.

    Returns:
        tuple of lists containing the reference names and reference lengths in
        corresponding order
    """
    with AlignmentFile(sam_path, ignore_truncation=True) as sam:
        contigs = list(sam.references)
        sizes = list(sam.lengths)
        indices = range(len(contigs))
        indices = sorted(indices, key=lambda idx: contigs[idx])
        return (
            [contigs[idx] for idx in indices],
            [sizes[idx] for idx in indices],
        )
