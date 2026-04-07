from reditools.alignment_manager import AlignmentManager


def setup_alignment_manager(
    file_list,
    min_read_quality,
    min_read_length,
    exclusions_file,
):
    """
    Create an AlignmentManager for REDItools.

    Parameters:
        options (namespace): Commandline arguments
        file_list (list): BAM file paths
        min_read_quality (int): Filter out reads with a MAPQ below threshold
        min_read_length (int): Filter out reads with length below threshold
        exclusions_file (str): Path to text file with read names to exclude

    Returns:
        AlignmentManager
    """
    sam_manager = AlignmentManager(
        ignore_truncation=True,
    )
    sam_manager.min_quality = min_read_quality
    sam_manager.min_length = min_read_length
    for sam in file_list:
        sam_manager.add_file(
            sam,
            exclusions_file,
        )
    return sam_manager
