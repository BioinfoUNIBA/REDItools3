from reditools.logger import Logger


def check_variants(options, bases):
    """
    Check whether specified edits are present.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """

    for variant in bases.variants:
        if variant in options.variants:
            return None
    return (
        'DISCARD COLUMN Requested edits {} not found',
        options.variants,
    )
