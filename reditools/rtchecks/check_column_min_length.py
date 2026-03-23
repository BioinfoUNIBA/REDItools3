from reditools.logger import Logger

def check_column_min_length(rtools, bases):
    """
    Check read depth.

    Parameters:
        rtools (REDItools): Object performing analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if the read depth is sufficient
    """
    if len(bases) < rtools.min_column_length:
        rtools.log(
            Logger.debug_level,
            'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
            len(bases),
            rtools.min_column_length,
        )
        return False
    return True
