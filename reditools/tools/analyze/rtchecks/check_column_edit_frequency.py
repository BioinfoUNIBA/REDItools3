"""Check if a position has a minimum number of total edits."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse

    from reditools.compiled_position import RTResult


class CheckColumnEditFrequency:
    """Check if a position has a minimum number of alternate bases.

    An alternate base is when there is at least one read representing a base
    other than the reference. For example, if the reference base is A, then
    the alternate possible bases are C, G, and T. If there are no detected
    edits, then the column edit frequency is 0. If there is one or more reads
    with a C at this position, than the column edit frequency would be 1. The
    maxmimum possible column edit frequency is 3.

    Attributes
    ----------
    min_edits : int
        The minimum required total edits.
    """

    def __init__(self, options: argparse.Namespace) -> None:
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

    def run_check(self, rtresult: RTResult) -> tuple | None:
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
        variant_no = len(rtresult.variants)
        if variant_no < self.min_edits:
            return (
                "DISCARDING COLUMN edits={} < {}",
                variant_no,
                self.min_edits,
            )
        return None
