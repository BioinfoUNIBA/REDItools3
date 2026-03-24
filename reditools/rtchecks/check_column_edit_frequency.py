from reditools.logger import Logger


def check_column_edit_frequency(rtools, bases):
    """
    Check the number of edits at the site.

    Parameters:
        rtools (REDItools): Object performing analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there are sufficient edits.
    """
    edits_no = len(bases) - bases['REF']
    if edits_no < rtools.min_edits:
        rtools.log(
            Logger.debug_level,
            'DISCARDING COLUMN edits={} < {}',
            edits_no,
            rtools.min_edits,
        )
        return False
    return True
