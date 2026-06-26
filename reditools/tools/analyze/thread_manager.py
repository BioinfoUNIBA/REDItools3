from __future__ import annotations

import argparse
import warnings
from multiprocessing import Process, Queue

from reditools.region import Region
from reditools.tools.analyze.redi_thread import redi_thread


class ThreadManager:
    def __init__(self, threads: int):
        self.nthreads = threads
        self.queue: Queue[tuple[Region, str] | None] = Queue()

    def start_threads(self, options: argparse.Namespace) -> None:
        self.processes = []
        for _ in range(self.nthreads):
            thread = Process(
                target=redi_thread,
                args=(options, self.queue),
            )
            thread.start()
            self.processes.append(thread)

    def kill_all_threads(self) -> None:
        for proc in self.processes:
            proc.kill()

    def await_finish(self) -> None:
        """Monitor progress of parallel analysis processes.

        Parameters
        ----------
        processes : list[Process]
            The list of worker processes.
        """
        is_running = True
        while is_running:
            is_running = False
            for process in self.processes:
                if process.exitcode == 1:
                    self.kill_all_threads()
                    raise Exception('Proccess died unexpectedly.')
                elif process.is_alive():
                    is_running = True

    def fill_queue(self, temp_filenames: list[tuple[Region, str]]) -> None:
        """
        Fill the input queue with genomic regions to be analyzed.

        Parameters
        ----------
        temp_filenames : list[tuple[Region, str]]
            List should contain tuples of Region obj representing range to
            analyze and a string representing the file to save results to.
        """
        if len(temp_filenames) < self.nthreads:
            warnings.warn(
                (
                    f"You have assigned more threads ({self.nthreads}) "
                    f"than there are genomic ranges ({len(temp_filenames)})."
                ),
                Warning,
            )

        for arg_tuple in temp_filenames:
            self.queue.put(arg_tuple)
        for _ in range(self.nthreads):
            self.queue.put(None)
