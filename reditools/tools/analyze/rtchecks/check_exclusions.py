from reditools.logger import Logger


def check_exclusions(rtools, bases):
    """
    Check if the bases object is in an excluded position.

    Parameters:
        namespace (namespace): Object running the analysis
        bases (CompiledPosition): Data for analysis

    Returns:
        (bool): True if the position is not excluded
    """
    in_exclusions = rtools.exclude_regions.contains(
        bases.contig,
        bases.position,
    )
    if in_exclusions:
        return ('DISCARD COLUMN in excluded region',)
