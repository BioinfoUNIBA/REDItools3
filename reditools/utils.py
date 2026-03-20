"""Miscellaneous utility functions."""

import csv
import os
import socket

from reditools.file_utils import open_stream


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
