"""Commandline tool for REDItools."""
import argparse
import sys
import traceback
from multiprocessing import Queue

from reditools.alignment_manager import AlignmentManager
from reditools.reditools import REDItools
from reditools.region import Region
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.setup_alignment_manager import \
    setup_alignment_manager
from reditools.tools.analyze.setup_rtools import setup_rtools
from reditools.tools.analyze.write_results import write_results


def analyze(
        options: argparse.Namespace,
        rtools: REDItools,
        sam_manager: AlignmentManager,
        region: Region,
        rtqc: RTChecks,
) -> str:
    rtresults = rtools.analyze(sam_manager, region)
    return write_results(
        rtresults,
        options.output_format,
        options.temp_dir,
        rtqc,
        rtools.log,
    )

def redi_thread(
        options: argparse.Namespace,
        in_queue: Queue,
        out_queue: Queue,
) -> bool:
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
    rtqc = RTChecks(options)
    while True:
        args = in_queue.get()
        if args is None:
            return True
        idx, region = args
        try:  # noqa: WPS229
            out_queue.put((
                idx,
                analyze(options, rtools, sam_manager, region, rtqc),
            ))
        except Exception as exc:
            if options.debug:
                traceback.print_exception(*sys.exc_info())
            sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
            sys.exit(1)
