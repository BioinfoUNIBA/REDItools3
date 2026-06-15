from reditools import file_utils
from reditools.alignment_manager import AlignmentManager


def setup_alignment_manager(
    file_list: list[str],
    min_read_quality: int,
    min_read_length: int,
    exclusions_file: str | None,
) -> AlignmentManager:
    """Initialize and configure an AlignmentManager object.

    Parameters
    ----------
    file_list : list[str]
        List of paths to alignment files (BAM/SAM).
    min_read_quality : int
        Minimum mapping quality for a read to be considered.
    min_read_length : int
        Minimum length for a read to be considered.
    exclusions_file : str, optional
        Path to a file containing read names to be excluded.

    Returns
    -------
    AlignmentManager
        A configured AlignmentManager instance.
    """
    if exclusions_file:
        exclude_set = set(file_utils.load_text_file(exclusions_file))
    else:
        exclude_set = None

    sam_manager = AlignmentManager(
        min_quality = min_read_quality,
        min_length = min_read_length,
        excluded_read_names=exclude_set,
    )
    for sam in file_list:
        sam_manager.add_file(
            sam,
        )
    return sam_manager
