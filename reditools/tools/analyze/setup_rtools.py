from reditools import file_utils, reditools
from reditools.logger import Logger


def setup_rtools(options):
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

    if options.variants:
        rtools.specific_edits = [_.upper() for _ in options.variants]

    if options.bed_file:
        for fname in options.bed_file:
            regions = file_utils.read_bed_file(fname)
            rtools.add_target_regions(regions)
    if options.exclude_regions:
        for fname in options.exclude_regions:
            regions = file_utils.read_bed_file(fname)
            rtools.add_exclude_regions(regions)
    if options.reference:
        rtools.add_reference(options.reference)

    if options.splicing_file:
        rtools.splice_positions = file_utils.load_splicing_file(
            options.splicing_file,
            options.splicing_span,
        )
        rtools.add_exclude_regions(regions)

    rtools.min_base_position = options.min_base_position
    rtools.max_base_position = options.max_base_position
    rtools.min_base_quality = options.min_base_quality

    rtools.min_column_length = options.min_read_depth
    rtools.min_edits = options.min_edits
    rtools.min_edits_per_nucleotide = options.min_edits_per_nucleotide
    rtools.max_alts = options.max_editing_nucleotides

    rtools.strand = options.strand
    rtools.strand_confidence_threshold = options.strand_confidence_threshold

    if options.strand_correction:
        rtools.use_strand_correction()

    return rtools
