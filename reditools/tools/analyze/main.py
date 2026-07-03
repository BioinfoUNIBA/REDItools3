from __future__ import annotations

import argparse
import sys
from functools import partial
from multiprocessing.context import TimeoutError
from multiprocessing.pool import Pool

from reditools.logger import Logger
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import parse_args
from reditools.tools.analyze.redi_thread import REDIThreadManager
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

def pool_error(pool: Pool, debug: bool, exc: Exception) -> None:
    """
    Terminates a multiprocessing Pool.

    Parameters
    ----------
    pool : Pool
        mutliprocessing Pool to terminate.
    debug : bool
        If True, raises the exception passed in the third argument.
    exc : Exception
        Exception responsible for the pool to terminate.
    """
    pool.terminate()
    if debug:
        raise exc.__cause__  # type: ignore[misc]
    sys.stderr.write(f'[ERROR] ({type(exc)}) {exc}\n')

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

    if options.threads > len(regions):
        sys.stderr.write(
            f"[WARNING] You have assigned {options.threads} threads, "
            f"But there are only {len(regions)} genomic range(s). "
            "Consider change the value of --window\n"
        )
        options.threads = len(regions)
    try:
        with Pool(
            options.threads,
            REDIThreadManager.init_thread,
            (options,),
        ) as pool:
            imap_iter = [
                pool.apply_async(
                    REDIThreadManager.analyze,
                    args=(region,),
                    error_callback=partial(pool_error, pool, options.debug),
                ) for region in regions
            ]
            pool.close()
            pool.join()
            temp_files = [_.get(1) for _ in imap_iter]
    except (TimeoutError, IndexError):
        sys.exit(1)

    concat_output(
        temp_files,
        options.output_file,
        'a' if options.append_file else 'w',
        options.encoding,
    )

    logger.log(Logger.info_level, 'Analyze Complete!')
