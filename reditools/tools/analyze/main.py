from __future__ import annotations

import argparse
import sys
import traceback
from multiprocessing import Pool

from reditools.logger import Logger
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import parse_args
from reditools.tools.analyze.redi_thread import REDIThread
from reditools.tools.analyze.region_args import region_args


def options_to_string(options: argparse.Namespace) -> str:
    """
    Convert argparse options to a comma-separated string of key:value pairs.

    Parameters
    ----------
    options : argparse.Namespace
        The parsed command line options.

    Returns
    -------
    str
        A string representation of the options.
    """
    return ", ".join(
        [f"{_}:{getattr(options, _)}" for _ in vars(options)],  # noqa: WPS421
    )

def setup_logger(options: argparse.Namespace) -> Logger:
    """
    Configure a logger based on the command line options.

    Parameters
    ----------
    options : argparse.Namespace
        The parsed command line options.

    Returns
    -------
    Logger
        The configured Logger object.
    """
    if options.debug:
        return Logger(Logger.debug_level)
    if options.verbose:
        return Logger(Logger.info_level)
    return Logger(Logger.silent_level)

def main() -> None:
    """
    The main entry point for the REDItools analyze command.
    """
    options = parse_args()

    logger = setup_logger(options)

    logger.log(logger.info_level, 'Starting REDItools')
    logger.log(
        logger.info_level,
        "Summary of command line parameters: {}",
        options_to_string(options),
    )

    options.encoding = 'utf-8'

    regions = region_args(options)
    # Re-implement thread count warning here
    try:
        with Pool(options.threads, REDIThread.init, (options,)) as pool:
            imap_iter = pool.imap(REDIThread.analyze, regions, 1)
            temp_files = [imap_iter.next() for _ in range(len(regions))]
    except Exception as exc:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')
        sys.exit(1)

    concat_output(
        temp_files,
        options.output_file,
        'a' if options.append_file else 'w',
        options.encoding,
    )

    logger.log(Logger.info_level, 'Analyze Complete!')
