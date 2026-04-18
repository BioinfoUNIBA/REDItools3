import argparse

from reditools import reditools
from reditools.logger import Logger


def setup_rtools(options: argparse.Namespace) -> reditools.REDItools:
    """
    Create a REDItools object.

    Parameters:
        options (namespace): Commandline arguments from argparse

    Returns:
        A configured REDItools object
    """
    rtools = reditools.REDItools()

    if options.debug:
        rtools.log_level = Logger.debug_level
    elif options.verbose:
        rtools.log_level = Logger.info_level

    if options.reference:
        rtools.add_reference(options.reference)

    rtools.min_base_position = options.min_base_position
    rtools.max_base_position = options.max_base_position
    rtools.min_base_quality = options.min_base_quality

    rtools.strand = options.strand
    rtools.strand_confidence_threshold = options.strand_confidence_threshold

    if options.strand_correction:
        rtools.use_strand_correction()

    return rtools
