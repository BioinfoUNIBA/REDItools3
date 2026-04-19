from __future__ import annotations

import sys
from multiprocessing import Process, Queue
from queue import Empty as EmptyQueueException


def check_dead(processes: list[Process]):
    """Check if any of the processes have failed.

    If a process has exited with code 1, all other processes are killed
    and the program exits.

    Parameters
    ----------
    processes : list[Process]
        The list of processes to monitor.
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
    """Monitor progress of parallel analysis processes.

    Parameters
    ----------
    processes : list[Process]
        The list of worker processes.
    out_queue : Queue[tuple[int, str]]
        The queue containing analysis result filenames and their indices.
    chunks : int
        The total number of work chunks.

    Returns
    -------
    list[str]
        A list of filenames containing analysis results, ordered by chunk index.
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
