from reditools.logger import Logger


def check_max_editing_nucleotides(options, bases):
    """
    Check that there are no more than a max number of alts.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    alts = bases.alts
    if len(alts) > options.max_editing_nucleotides:
        return (
            'DISCARD COLUMN alts={} > {}',
            len(alts),
            options.max_editing_nucleotides,
        )
