from reditools.logger import Logger


def check_target_positions(rtools, bases):
    """
    Check if the bases object is in a target region.

    Parameters:
        namespace (namespace): Object running the analysis
        bases (CompiledPosition): Data for analysis

    Returns:
        (bool): True if the position is in a target region
    """
    if not rtools.target_regions.contains(
            bases.contig,
            bases.position,
    ):
        return ('DISCARD COLUMN not in target regions',)
