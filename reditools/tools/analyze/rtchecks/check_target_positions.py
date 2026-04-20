import argparse

from reditools import file_utils
from reditools.compiled_position import RTResult
from reditools.region_collection import RegionCollection


class CheckTargetPositions:
    """Check if a position is within target regions.

    Attributes
    ----------
    regions : RegionCollection
        The collection of target regions.
    """

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckTargetPositions.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing target BED files.
        """
        self.regions = RegionCollection()
        self.regions.add_regions(file_utils.read_bed_file(*options.bed_file))

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        """Check if this check is required based on options.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.

        Returns
        -------
        bool
            True if bed_file option is provided, False otherwise.
        """
        return options.bed_file is not None

    def run_check(self, bases: RTResult) -> None | tuple:
        """Run the check on a specific position.

        Parameters
        ----------
        bases : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        None | tuple
            None if the position is within target regions, a tuple with
            error message otherwise.
        """
        if not self.regions.contains(bases.contig, bases.position):
            return ('DISCARD COLUMN not in target regions',)
        return None
