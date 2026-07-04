import argparse
import sys
import traceback
from functools import partial
from multiprocessing.context import TimeoutError
from multiprocessing.pool import Pool

from reditools.tools.analyze.redi_thread import REDIThreadManager
from reditools.tools.analyze.temp_file_manager import TempFileManager


def run_pool(
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
                    error_callback=partial(terminate_pool, pool, options.debug),
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

def terminate_pool(
    pool: Pool,
    debug: bool,
    exc: Exception,
) -> None:
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


