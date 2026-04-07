"""Commandline tool for REDItools."""

from queue import Empty as EmptyQueueException
from .check_dead import check_dead


def monitor(processes, out_queue, chunks):
    """
    Monitor parallel REDItools jobs.

    Parameters:
        processes (list): Threads
        out_queue (Queue): Output of threads
        chunks (int): Number of chunks for analysis

    Returns:
        list: Temporary files containing the output of each chunk.
    """
    tfs = [None for _ in range(chunks - len(processes))]

    for prc in processes:
        prc.start()

    while None in tfs:
        try:
            idx, fname = out_queue.get(block=False, timeout=1)
            tfs[idx] = fname
        except EmptyQueueException:
            check_dead(processes)
    return tfs
