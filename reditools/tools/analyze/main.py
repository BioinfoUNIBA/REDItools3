from __future__ import annotations

import os
import sys
import tempfile
import traceback
from functools import partial
from multiprocessing.context import TimeoutError
from multiprocessing.pool import Pool
from typing import TYPE_CHECKING

from reditools.logger import Logger
from reditools.tools.analyze.parse_args import json_args, parse_args
from reditools.tools.analyze.redi_thread import REDIThreadManager
from reditools.tools.analyze.region_args import region_args
from reditools.tools.analyze.temp_file_manager import TempFileManager
from reditools import file_utils

if TYPE_CHECKING:
    import argparse

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

def analyze(
    options: argparse.Namespace,
    temp_filemanager: TempFileManager,
) -> bool:
    """
    Create a pool of threads and analyze the data.

    Parameters
    ----------
    options : argparse.Namespace
        CLI arguments.
    temp_filemanager : TempFileManager
        Regions to analyze and files to save to.

    Returns
    -------
    bool
        True if the analysis completes successfully, False otherwise.
    """
    try:
        with Pool(
            options.threads,
            REDIThreadManager.init_thread,
            (options,),
        ) as pool:
            imap_iter = [
                pool.apply_async(
                    REDIThreadManager.analyze,
                    args=(region, filename),
                    error_callback=partial(pool_error, pool, options.debug),
                ) for region, filename in temp_filemanager
            ]
            pool.close()
            pool.join()
            [_.get(1) for _ in imap_iter]
    except TimeoutError:
        return False
    except Exception:
        if options.debug:
            traceback.print_exception(*sys.exc_info())
        return False
    return True

def main() -> None:
    """
    The main entry point for the REDItools analyze command.
    """
    options = parse_args.parse_args()

    logger = setup_logger(options)

    if options.resume:
        logger.log(
            logger.info_level,
            (
                'Resuming REDItools from directory "{}". Using parameters '
                'from previous run. All other command line options will be '
                'ignored.'
            ),
            options.temp_dir,
        )
        temp_dir = options.temp_dir
    else:
        logger.log(logger.info_level, 'Starting REDItools')
        temp_dir = file_utils.make_temp_dir(
            prefix='reditools_',
            dir=options.temp_dir,
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

    if options.resume:
        temp_file_manager = TempFileManager(temp_dir)
    else:
        temp_file_manager = TempFileManager(temp_dir, region_args(options))

    if options.threads > len(temp_file_manager):
        sys.stderr.write(
            f"[WARNING] You have assigned {options.threads} threads, "
            f"But there are only {len(temp_file_manager)} genomic range(s). "
            "Consider change the value of --window\n"
        )
        options.threads = len(temp_file_manager)

    if not analyze(options, temp_file_manager):
        sys.exit(1)

    temp_file_manager.concat(
        options.output_file,
        'a' if options.append_file else 'w',
    )
    temp_file_manager.cleanup()

    logger.log(Logger.info_level, 'Analyze Complete!')
