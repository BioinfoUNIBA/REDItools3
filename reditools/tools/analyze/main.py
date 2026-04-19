"""Commandline tool for REDItools."""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Process, Queue

from reditools.logger import Logger
from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.monitor import monitor
from reditools.tools.analyze.parse_args import parse_args
from reditools.tools.analyze.redi_thread import redi_thread
from reditools.tools.analyze.region_args import region_args


def options_to_string(options: argparse.Namespace) -> str:
    return ", ".join(
        [f"{_}:{getattr(options, _)}" for _ in vars(options)],  # noqa: WPS421
    )

def setup_logger(options: argparse.Namespace) -> Logger:
    if options.debug:
        return Logger(Logger.debug_level)
    if options.verbose:
        return Logger(Logger.info_level)
    return Logger(Logger.silent_level)

def fill_queue(options: argparse.Namespace) -> Queue[tuple[int, Region] | None]:
    in_queue: Queue[tuple[int, Region] | None] = Queue()
    try:
        for _ in enumerate(region_args(options)):  # noqa: WPS468
            in_queue.put(_)
    except FileNotFoundError as exc:
        sys.stderr.write(f'[ERROR] {exc}\n')
        sys.exit(1)

    # Check thread count
    if in_queue.qsize() < options.threads:
        sys.stderr.write(
            "[WARNING] You have assigned more threads "
            f"({options.threads}) than there are genomic ranges "
            f"({in_queue.qsize()})\n",
        )
        options.threads = in_queue.qsize()
    for _ in range(options.threads):
        in_queue.put(None)
    return in_queue

def main():
    """Perform RNA editing analysis."""
    options = parse_args()

    logger = setup_logger(options)

    logger.log(logger.info_level, 'Starting REDItools')
    logger.log(
        logger.info_level,
        "Summary of command line parameters: {}",
        options_to_string(options),
    )

    options.output_format = {'delimiter': '\t', 'lineterminator': '\n'}
    options.encoding = 'utf-8'

    in_queue = fill_queue(options)

    # Start parallel jobs
    out_queue = Queue()
    processes = []
    for _ in range(options.threads):
        processes.append(Process(
            target=redi_thread,
            args=(options, in_queue, out_queue),
        ))

    concat_output(
        monitor(processes, out_queue, in_queue.qsize()),
        options.output_file,
        'a' if options.append_file else 'w',
        options.encoding,
        **options.output_format)

    logger.log(Logger.info_level, 'Analyze Complete!')
