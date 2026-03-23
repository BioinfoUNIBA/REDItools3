from reditools.logger import Logger

def check_multiple_alts(rtools, bases):
    """
    Check that there is, at most, one alternate base.

    Parameters:
        rtools (REDItools): Object running the analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there is zero or one alt
    """
    alts = bases.get_variants()
    if len(alts) > 1:
        rtools.log(
            Logger.debug_level,
            'DISCARD COLUMN alts={} > 1',
            len(alts),
        )
        return False
    return True
