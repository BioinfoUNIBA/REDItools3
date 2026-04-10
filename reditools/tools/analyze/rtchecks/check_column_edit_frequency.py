from reditools.logger import Logger


def check_column_edit_frequency(options, bases):
    """
    Check the number of edits at the site.

    Parameters:
        options (namespace): Analyze tool options
        bases (CompiledPosition): Base position under analysis

    Returns:
        None if QC passed, else debug message (tuple)
    """
    edits_no = len(bases) - bases['REF']
    if edits_no < options.min_edits:
        return (
            'DISCARDING COLUMN edits={} < {}',
            edits_no,
            options.min_edits,
        )
