from reditools.logger import Logger

__all__ = ('check_column_min_edits')
_bases = ('A', 'T', 'C', 'G')


def check_column_min_edits(rtools, bases):
    """
    Check that there are sufficient edit events for each base.

    Parameters:
        rtools (REDItools): Object performing analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there are sufficient edits
    """
    for base in _bases:
        if base != bases.ref and \
                rtools.min_edits_per_nucleotide > bases[base] > 0:
            rtools.log(
                Logger.debug_level,
                'DISCARDING COLUMN edits={} < {}',
                bases[base],
                rtools.min_edits_per_nucleotide,
            )
            return False
    return True
