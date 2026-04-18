"""Commandline tool for REDItools."""
from __future__ import annotations

import sys
from multiprocessing import Process, Queue
from queue import Empty as EmptyQueueException


def check_dead(processes: list[Process]):
    """
    Look through processes to determine if any have died unexpectedly.

    If any process has an exit code of 1, this method will terminate all other
    processes and then exit with code 1.

    Parameters:
        processes (list): Processes to check
    """
    for proc in processes:
        if proc.exitcode == 1:
            for to_kill in processes:
                to_kill.kill()
            sys.stderr.write('[ERROR] Killing job\n')
            sys.exit(1)

def monitor(
        processes: list[Process],
        out_queue: Queue[tuple[int, str]],
        chunks: int,
) -> list[str]:
    """
    Monitor parallel REDItools jobs.

    Parameters:
        processes (list): Threads
        out_queue (Queue): Output of threads
        chunks (int): Number of chunks for analysis

    Returns:
        list: Temporary files containing the output of each chunk.
    """
    tfs = ['' for _ in range(chunks - len(processes))]

    for prc in processes:
        prc.start()

    while '' in tfs:
        try:
            idx, fname = out_queue.get(block=False, timeout=1)
        except EmptyQueueException:
            check_dead(processes)
        else:
            tfs[idx] = fname
    return tfs
