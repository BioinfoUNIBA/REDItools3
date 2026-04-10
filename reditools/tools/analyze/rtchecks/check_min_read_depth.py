from reditools.logger import Logger


def check_min_read_depth(options, bases):
    """
    Checks whether there is sufficient read coverage.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    if len(bases) < options.min_read_depth:
        return (
            'DISCARDING COLUMN {} [MIN_READ_DEPTH={}]',
            len(bases),
            options.min_read_depth,
        )
