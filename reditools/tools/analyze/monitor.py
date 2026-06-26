from __future__ import annotations

import sys
from multiprocessing import Process


def kill_all(processes: list[Process]) -> None:
    for proc in processes:
        proc.kill()

def check_dead(processes: list[Process]) -> None:
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
) -> None:
    """Monitor progress of parallel analysis processes.

    Parameters
    ----------
    processes : list[Process]
        The list of worker processes.
    """
    for prc in processes:
        prc.start()

    is_running = True
    while is_running:
        is_running = False
        for proc in processes:
            if proc.exitcode == 1:
                kill_all(processes)
                sys.stderr.write('[ERROR] Killing job\n')
                sys.exit(1)
            elif proc.is_alive():
                is_running = True
