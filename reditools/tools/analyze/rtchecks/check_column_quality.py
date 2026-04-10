from reditools.logger import Logger


# Really shouldn't use this one. I have to compute mean_q anyway
def check_column_quality(options, bases):
    """
    Check mean quality of the position.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    if bases:
        mean_q = sum(bases.qualities) / len(bases)
    else:
        mean_q = 0
    if mean_q < options.min_read_quality:
        return (
            Logger.debug_level,
            'DISCARD COLUMN mean_quality={} < {}',
            mean_q,
            options.min_read_quality,
        )
