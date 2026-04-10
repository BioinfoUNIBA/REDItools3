from reditools.logger import Logger


def check_target_positions(options, bases):
    """
    Check if the bases object is in a target region.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    if not options.target_regions.contains(
            bases.contig,
            bases.position,
    ):
        return ('DISCARD COLUMN not in target regions',)
