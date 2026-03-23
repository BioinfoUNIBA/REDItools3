from reditools.logger import Logger

def check_target_positions(rtools, bases):
    """
    Check if the bases object is in a target region.

    Parameters:
        rtools (REDItools): Object running the analysis
        bases (CompiledPosition): Data for analysis

    Returns:
        (bool): True if the position is in a target region
    """
    in_targets = rtools.target_regions.contains(
        bases.contig,
        bases.position,
    )
    if not in_targets:
        rtools.log(
            Logger.debug_level,
            'DISCARD COLUMN not in target regions',
        )
    return in_targets
