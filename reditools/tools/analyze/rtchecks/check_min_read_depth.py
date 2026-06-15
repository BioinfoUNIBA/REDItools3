import argparse

from reditools.compiled_position import RTResult


class CheckMinReadDepth:
    """Check if a position has minimum read depth.

    Attributes
    ----------
    min_read_depth : int
        The minimum required read depth.
    """

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckMinReadDepth.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing min_read_depth.
        """
        self.min_read_depth = options.min_read_depth

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
            True if min_read_depth > 1, False otherwise.
        """
        return options.min_read_depth > 1

    def run_check(self, rtresult: RTResult) -> None | tuple:
        """Run the check on a specific position.

        Parameters
        ----------
        rtresult : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        None | tuple
            None if read depth is sufficient, a tuple with error message
            otherwise.
        """
        if len(rtresult) < self.min_read_depth:
            return (
                'DISCARDING COLUMN {} [MIN_READ_DEPTH={}]',
                len(rtresult),
                self.min_read_depth,
            )
        return None
