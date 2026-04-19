import argparse

from reditools import file_utils
from reditools.compiled_position import RTResult
from reditools.region_collection import RegionCollection


class CheckExclusions:
    """Check if a position is within excluded regions.

    Attributes
    ----------
    regions : RegionCollection
        The collection of excluded regions.
    """

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckExclusions.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing excluded regions BED files.
        """
        self.regions = RegionCollection()
        self.regions.add_regions(file_utils.read_bed_file(*options.exclude_regions))

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
            True if exclude_regions option is provided, False otherwise.
        """
        return options.exclude_regions is not None

    def run_check(self, bases: RTResult) -> None | tuple:
        """Run the check on a specific position.

        Parameters
        ----------
        bases : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        None | tuple
            None if the position is not in excluded regions, a tuple with
            error message otherwise.
        """
        if self.regions.contains(bases.contig, bases.position):
            return ('DISCARD COLUMN in excluded region',)
        return None
