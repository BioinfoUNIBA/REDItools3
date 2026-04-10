from reditools.logger import Logger


def check_exclusions(options, bases):
    """
    Check if the bases object is in an excluded position.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    in_exclusions = options.exclude_regions.contains(
        bases.contig,
        bases.position,
    )
    if in_exclusions:
        return ('DISCARD COLUMN in excluded region',)
