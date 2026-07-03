from __future__ import annotations

import json
import os
import sys
import tempfile
import traceback
from functools import partial
from multiprocessing.context import TimeoutError
from multiprocessing.pool import Pool
from typing import TYPE_CHECKING

from reditools.logger import Logger
from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import json_args, parse_args
from reditools.tools.analyze.redi_thread import REDIThreadManager
from reditools.tools.analyze.region_args import region_args

if TYPE_CHECKING:
    import argparse

json_windows_file = 'tempfile_map.json'

def make_temp_dir(prefix: str | None=None, dir: str | None=None) -> str:
    """
    Creates a folder.

    Parameters
    ----------
    prefix : str
        Filename prefix.
    dir : str
        Path to folder parent.

    Returns
    -------
    str
        Path to the folder.
    """
    with tempfile.NamedTemporaryFile(prefix=prefix, dir=dir) as stream:
        valid_name = stream.name
    os.mkdir(valid_name)
    return valid_name

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

def get_temp_filenames_list(
    options: argparse.Namespace,
    temp_dir: str,
) -> list[tuple[Region, str]]:
    """
    Creates a list of genomic ranges and filenames to store analysis results.

    Parameters
    ----------
    options : argparse.Namespace
        CLI options.
    temp_dir : str
        Folder to store analysis result files.

    Returns
    -------
    list[tuple[Region, str]]
        Each list element will contain a tuple of the genomic range of the
        analysis segment and the file path that will eventually contain
        the analysis results.
    """
    if options.resume:
        with open(os.path.join(temp_dir, json_windows_file), 'r') as stream:
            temp_filenames = [
                (Region.from_string(region), filename)
                for region, filename in json.load(stream)
            ]
    else:
        temp_filenames = [
            (
                region,
                tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name,
            )
            for region in region_args(options)
        ]
        with open(os.path.join(temp_dir, json_windows_file), 'w') as stream:
            json.dump(temp_filenames, stream, default=str)
    return temp_filenames

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

def cleanup_temp_files(
    temp_dir: str,
    temp_filenames: list[tuple[Region, str]],
) -> None:
    """
    Deletes the temporary files and directory made by the tool.

    Parameters
    ----------
    temp_dir : str
        Path to analyze temporary directory.
    temp_filenames : list[tuple[Region, str]]
        List of temporary analysis files.
    """
    for _, temp_file in temp_filenames:
        os.remove(f'{temp_file}.done')

    for temp_file in (json_args.json_args_filename, json_windows_file):
        os.remove(os.path.join(temp_dir, temp_file))
    try:
        os.rmdir(temp_dir)
    except OSError as exc:
        sys.stderr.write(
            f'[WARNING] Could not delete temporary files directory {temp_dir}. '
            f'{exc}\n'
        )

def analyze(
    options: argparse.Namespace,
    temp_filenames: list[tuple[Region, str]],
) -> bool:
    """
    Create a pool of threads and analyze the data.

    Parameters
    ----------
    options : argparse.Namespace
        CLI arguments.
    temp_filenames : list[tuple[Region, str]]
        Regions to analyze and the files to store them in.

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
                ) for region, filename in temp_filenames
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
        temp_dir = make_temp_dir(
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

    options.encoding = 'utf-8'

    temp_filenames = get_temp_filenames_list(options, temp_dir)

    if options.threads > len(temp_filenames):
        sys.stderr.write(
            f"[WARNING] You have assigned {options.threads} threads, "
            f"But there are only {len(temp_filenames)} genomic range(s). "
            "Consider change the value of --window\n"
        )
        options.threads = len(temp_filenames)

    if not analyze(options, temp_filenames):
        sys.exit(1)

    concat_output(
        [_[1] for _ in temp_filenames],
        options.output_file,
        'a' if options.append_file else 'w',
        options.encoding,
    )

    cleanup_temp_files(temp_dir, temp_filenames)

    logger.log(Logger.info_level, 'Analyze Complete!')
