from reditools.logger import Logger


def check_variants(namespace, bases):
    """
    Check whether specified edits are present.

    Parameters:
        namespace (namespace): Object running the analysis
        bases (CompiledPosition): Base position under analysis

    Returns:
        (bool): True if there specified edits are detected.
    """

    for variant in bases.variants:
        if variant in namespace.variants:
            return None
    return (
        'DISCARD COLUMN Requested edits {} not found',
        namespace.variants,
    )
