"""Commandline tool for REDItools."""
import traceback
import sys

from reditools.tools.analyze.setup_rtools import setup_rtools
from reditools.tools.analyze.setup_alignment_manager import setup_alignment_manager
from reditools.tools.analyze.write_results import write_results


def run_proc(options, in_queue, out_queue):
    """
    Analyze a genomic segment using REDItools.

    Parameters:
        options (namesapce): Configuration options from argparse for REDItools
        in_queue (Queue): Queue of input arguments for analysis
        out_queue (Queue): Queue to store paths to analysis results

    Returns:
        bool: True if the in_queue is empty
    """
    rtools = setup_rtools(options)
    sam_manager = setup_alignment_manager(
        options.file,
        options.min_read_quality,
        options.min_read_length,
        options.exclude_reads,
    )
    while True:
        args = in_queue.get()
        if args is None:
            return True
        idx, region = args
        try:  # noqa: WPS229
            rtresults = rtools.analyze(sam_manager, region)
            out_queue.put((idx, write_results(
                rtresults,
                options.file,
                options.output_format,
                options.temp_dir,
            )))
        except Exception as exc:
            if options.debug:
                traceback.print_exception(*sys.exc_info())
            sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
            sys.exit(1)

