from reditools.logger import Logger

_bases = ('A', 'T', 'C', 'G')


def check_column_min_edits(options, bases):
    """
    Check that there are sufficient edit events for each base.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    for base in _bases:
        if base != bases.reference and \
                0 < bases[base] < options.min_edits_per_nucleotide:
            return (
                'DISCARDING COLUMN edits={} < {}',
                bases[base],
                options.min_edits_per_nucleotide,
            )
