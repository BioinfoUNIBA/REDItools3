#!/usr/bin/env python

'''
Miscellaneous utility functions
'''

import os
import socket
import csv
from gzip import open as gzip_open
from collections import defaultdict
from sortedcontainers import SortedSet

def open_stream(path, mode="rt", encoding="utf-8"):
    '''
    Opens a input or output stream from a file, accounting for gzip.

    Parameters:
        path (str): Path to file for reading or writing
        mode (str): File mode

        gzip (bool): Whether the file is or should be gzipped

    Returns:
        TextIOWrapper to the file
    '''
    if path.endswith("gz"):
        return gzip_open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)

def read_bed_file(path):
    '''
    Returns an iterator for a BED file.

    Parameters:
        path (str): Path to a BED file for reading.

    Returns:
        Iterator of BED file contents.
    '''
    stream = open_stream(path)
    return csv.reader(stream, delimiter="\t")

def enumerate_positions(regions):
    '''
    Converts a list of regions into a list of individual positions.

    Parameters:
        regions (list): A list of iterables. Each element must start
                        with a contig and start position. End position
                        is optional. Additional values will be ignored.

    Returns:
        Dictionary of contigs to SortedSets enumerating the individual positions.
    '''
    positions = defaultdict(SortedSet)
    for reg in regions:
        contig = reg[0]
        start = int(reg[1]) - 1
        end = start + 1 if len(reg) < 3 else int(reg[2])
        positions[contig] |= range(start, end)
    return positions

def get_hostname_string():
    '''
    Retrieves the machine hostname, ip, and proccess ID.

    Returns:
        String in the format "hostname|ip|pid"
    '''
    hostname = socket.gethostname()
    ip_addr = socket.gethostbyname(hostname)
    pid = os.getpid()
    return f"{hostname}|{ip_addr}|{pid}"

def check_list(functions, **kwargs):
    '''
    Runs through a list of functions, determining if any return False.

    Parameters:
        functions (list): A list of function references
        **kwargs: Any arguments to be passed to the members of functions

    Returns:
        False if any function in check_list returns False, else True
    '''
    for check in functions:
        if not check(**kwargs):
            return False
    return True
