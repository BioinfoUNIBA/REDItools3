from reditools import file_utils
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

    if exclusions_file:
        with file_utils.open(exclusions_file, 'r') as stream:
            exclude_set = set(_.strip() for _ in stream)
    else:
        exclude_set = None

    sam_manager = AlignmentManager(
        ignore_truncation=True,
        exclude_set=exclude_set
    )
    sam_manager.min_quality = min_read_quality
    sam_manager.min_length = min_read_length
    for sam in file_list:
        sam_manager.add_file(
            sam,
        )
    return sam_manager
