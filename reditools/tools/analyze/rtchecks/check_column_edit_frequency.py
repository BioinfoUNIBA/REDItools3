import argparse

from reditools.compiled_position import RTResult


class CheckColumnEditFrequency:
    """Check if a position has a minimum number of total edits.

    Attributes
    ----------
    min_edits : int
        The minimum required total edits.
    """

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckColumnEditFrequency.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing min_edits.
        """
        self.min_edits = options.min_edits

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
            True if min_edits > 0, False otherwise.
        """
        return options.min_edits > 0

    def run_check(self, rtresult: RTResult) -> None | tuple:
        """Run the check on a specific position.

        Parameters
        ----------
        rtresult : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        None | tuple
            None if total edits are sufficient, a tuple with error message
            otherwise.
        """
        edits_no = len(rtresult) - rtresult['REF']
        if edits_no < self.min_edits:
            return (
                'DISCARDING COLUMN edits={} < {}',
                edits_no,
                self.min_edits,
            )
        return None
