from reditools.logger import Logger


def check_max_alts(rtools, bases):
    """
    Check that there are no more than a max number of alts.

    Parameters:
        rtools (REDItools): Object running the analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there are n or fewer alts
    """
    alts = bases.alts
    if len(alts) > rtools.max_alts:
        rtools.log(
            Logger.debug_level,
            'DISCARD COLUMN alts={} > {}',
            len(alts),
            rtools.max_alts,
        )
        return False
    return True
