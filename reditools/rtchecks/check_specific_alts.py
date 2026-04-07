from reditools.logger import Logger


def check_specific_alts(rtools, bases):
    """
    Check whether specified edits are present.

    Parameters:
        rtools (REDItools): Object running the analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there specified edits are detected.
    """

    for variant in bases.variants:
        if variant in rtools._specific_edits:
            return True
    rtools.log(
        Logger.debug_level,
        'DISCARD COLUMN Requested edits {} not found',
        rtools._specific_edits,
    )
    return False
