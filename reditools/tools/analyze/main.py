from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile

from reditools.logger import Logger
from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import json_args, parse_args
from reditools.tools.analyze.region_args import region_args
from reditools.tools.analyze.thread_manager import ThreadManager

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

def run_analysis(options: argparse.Namespace, temp_dir: str) -> None:
    """
    Starts the actual REDItools analysis.

    This method makes use of the ThreadManager class to handle
    the analysis process.

    Parameters
    ----------
    options : argparse.Namespace
        CLI arguments.
    temp_dir : str
        Folder to store temporary files.
    """
    temp_filenames = get_temp_filenames_list(options, temp_dir)

    thread_manager = ThreadManager(options.threads)
    thread_manager.fill_queue(temp_filenames)
    thread_manager.start_threads(options)
    thread_manager.await_finish()
    concat_output(
        [_[1] for _ in temp_filenames],
        options.output_file,
        'a' if options.append_file else 'w',
        'utf-8',
    )

    for _, temp_file in temp_filenames:
        os.remove(f'{temp_file}.done')

    for json_file in (json_args.json_args_filename, json_windows_file):
        os.remove(os.path.join(temp_dir, json_file))
    try:
        os.rmdir(temp_dir)
    except OSError as exc:
        sys.stderr.write(
            f'[WARNING] Could not delete temporary files directory {temp_dir}. '
            f'{exc}\n'
        )

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
   
    run_analysis(options, temp_dir) 


    logger.log(Logger.info_level, 'Analyze Complete!')
