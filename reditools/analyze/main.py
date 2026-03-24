"""Commandline tool for REDItools."""

import sys

from multiprocessing import Process, Queue
from reditools import file_utils
from reditools.region import Region

from .concat_output import concat_output
from .monitor import monitor
from .parse_options import parse_options
from .region_args import region_args
from .run_proc import run_proc


def main():
    """Perform RNA editing analysis."""
    options = parse_options()

    is_verbose = options.debug or options.verbose

    if is_verbose:
        options_string = ", ".join(
            [f"{_}:{getattr(options, _)}" for _ in vars(options)],
        )
        sys.stderr.write(
            "Starting REDItools\n"
            f"Summary of command line parameters: {options_string}\n",
        )

    options.output_format = {'delimiter': '\t', 'lineterminator': '\n'}
    options.encoding = 'utf-8'
    if options.exclude_reads:
        options.exclude_reads = file_utils.load_text_file(
            options.exclude_reads,
        )

    # Put analysis chunks into queue
    if options.region is None:
        region = None
    else:
        region = Region.from_string(options.region, options.file[0])

    try:
        regions = region_args(
            options.file[0],
            region,
            options.window,
        )
    except FileNotFoundError as e:
        sys.stderr.write(f'[ERROR] {e}\n')
        sys.exit(1)

    # Check thread count
    if len(regions) < options.threads:
        sys.stderr.write(
            "[WARNING] You have assigned more threads "
            f"({options.threads}) than there are genomic ranges "
            f"({len(regions)})\n",
        )
        options.threads = len(regions)

    in_queue = Queue()
    for args in enumerate(regions):
        in_queue.put(args)
    for _ in range(options.threads):
        in_queue.put(None)

    # Start parallel jobs
    out_queue = Queue()
    processes = [
        Process(
            target=run_proc,
            args=(options, in_queue, out_queue),
        ) for _ in range(options.threads)
    ]
    if is_verbose:
        sys.stderr.write(
            "All processes complete. Concatenating temporary files.\n",
        )

    concat_output(
        monitor(processes, out_queue, in_queue.qsize()),
        options.output_file,
        'a' if options.append_file else 'w',
        options.encoding,
        **options.output_format)

    if is_verbose:
        sys.stderr.write("Analaysis Complete!\n")
