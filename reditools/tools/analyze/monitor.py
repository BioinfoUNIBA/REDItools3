"""Commandline tool for REDItools."""
import sys
from queue import Empty as EmptyQueueException

def check_dead(processes):
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
