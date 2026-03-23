from reditools.logger import Logger

def check_is_none(rtools, bases):
    """
    Check if the bases object is None.

    Parameters:
        rtools (REDItools): Object running the analysis
        bases (CompiledPosition): Data for analysis

    Returns:
        (bool): True if bases is not None
    """
    if bases is None:
        rtools.log(
            Logger.debug_level,
            'DISCARD COLUMN no reads',
        )
        return False
    return True
