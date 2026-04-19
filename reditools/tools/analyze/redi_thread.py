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
    """Analyze a specific genomic region.

    Parameters
    ----------
    options : argparse.Namespace
        The command-line options.
    rtools : REDItools
        The REDItools analysis engine.
    sam_manager : AlignmentManager
        The alignment file manager.
    region : Region
        The genomic region to analyze.
    rtqc : RTChecks
        The quality control checks to apply.

    Returns
    -------
    str
        The path to the temporary file containing the results.
    """
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
    """Worker thread function for parallel REDItools analysis.

    Parameters
    ----------
    options : argparse.Namespace
        The command-line options.
    in_queue : Queue
        The queue containing genomic regions to analyze.
    out_queue : Queue
        The queue to put analysis results into.

    Returns
    -------
    bool
        True when the worker has finished processing all regions.
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
