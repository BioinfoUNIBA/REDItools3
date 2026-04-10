from reditools.logger import Logger


def check_max_editing_nucleotides(rtools, bases):
    """
    Check that there are no more than a max number of alts.

    Parameters:
        namespace (namespace): Object running the analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there are n or fewer alts
    """
    alts = bases.alts
    if len(alts) > rtools.max_editing_nucleotides:
        return (
            'DISCARD COLUMN alts={} > {}',
            len(alts),
            rtools.max_editing_nucleotides,
        )
