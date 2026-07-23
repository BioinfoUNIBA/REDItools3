"""REDItools analyze tool entry point form CLI."""
from __future__ import annotations

import resource
import sys
from typing import TYPE_CHECKING

from reditools import file_utils
from reditools.logger import Logger
from reditools.tools.analyze.parse_args import json_args, parse_args
from reditools.tools.analyze.redi_pool import run_pool
from reditools.tools.analyze.region_args import region_args
from reditools.tools.analyze.temp_file_manager import TempFileManager

if TYPE_CHECKING:
    import argparse

def main() -> None:
    """Begin REDItools analyze from CLI."""
    options = parse_args.parse_args()

    logger = setup_logger(options)

    if options.resume:
        logger.log(
            logger.info_level,
            (
                "Resuming REDItools from directory '{}'. Using parameters "
                "from previous run. All other command line options will be "
                "ignored."
            ),
            options.temp_dir,
        )
        temp_dir = options.temp_dir
    else:
        logger.log(logger.info_level, "Starting REDItools")
        temp_dir = file_utils.make_dir(
            prefix="reditools_",
            dirname=options.temp_dir,
        )
        json_args.args_to_json(options, temp_dir)

    logger.log(
        logger.info_level,
        "Summary of command line parameters: {}",
        parse_args.args_to_string(options),
    )

    logger.log(
        logger.info_level,
        "Temporary files will be written to {}",
        temp_dir,
    )

    if analyze(options, temp_dir):
        logger.log(Logger.info_level, "Analyze Complete!")
    else:
        sys.exit(1)

def setup_logger(options: argparse.Namespace) -> Logger:
    """Configure a logger based on the command line options.

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

def analyze(
    options: argparse.Namespace,
    temp_dir: str,
) -> bool:
    """Perform the REDItools analysis.

    Parameters
    ----------
    options : argparse.Namespace
        CLI arguments.
    temp_dir : str
        Where to save temporary files.

    Returns
    -------
    bool
        True if the analysis completed successfully, False otherwise.
    """
    with TempFileManager(
        temp_dir,
        None if options.resume else region_args(options),
    ) as temp_file_manager:
        options.threads = check_thread_limits(
            options.threads,
            len(temp_file_manager),
        )

        if not run_pool(options, temp_file_manager, options.dna_sample is not None):
            return False

        temp_file_manager.concat(
            options.output_file,
            "a" if options.append_file else "w",
        )
    return True

def check_thread_limits(n_threads: int, n_regions: int) -> int:
    """Limit the thread count based on user inputs and system limits.

    Thread count is reduced to the minimum of the total regions
    to analyze and the open file limit of the system.

    Parameters
    ----------
    n_threads : int
        Number of threads from CLI arguments
    n_regions : int
        Number of regions to analyze

    Returns
    -------
    int
        Maximum usable threads
    """
    n_file_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    if n_file_limit < n_regions and n_file_limit < n_threads:
        sys.stderr.write(
            f"[WARNING] Reducing your thread count to {n_file_limit} as your "
            f"system only allows a maximum of {n_file_limit} files open at "
            "a time. Consider increasing your system's open file limit.\n",
        )
        return n_file_limit
    if n_regions < n_threads:
        sys.stderr.write(
            f"[WARNING] You have assigned {n_threads} threads, "
            f"But there are only {n_regions} genomic range(s). "
            "Consider changing the value of --window\n",
        )
        return n_regions
    return n_threads
