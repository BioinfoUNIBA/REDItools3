import argparse

from reditools.compiled_position import RTResult


class CheckColumnMinEdits:
    """Check if a position has a minimum number of edits per nucleotide.
    Specifically, checks that all non-zero, non-reference bases pass a
    given threshold.

    Attributes
    ----------
    min_edits_per_nucleotide : int
        The minimum required edits per nucleotide.
    """

    _bases = ('A', 'T', 'C', 'G')

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckColumnMinEdits.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing min_edits_per_nucleotide.
        """
        self.min_edits_per_nucleotide = options.min_edits_per_nucleotide

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
            True if min_edits_per_nucleotide > 0, False otherwise.
        """
        return options.min_edits_per_nucleotide > 0

    def run_check(self, bases: RTResult) -> tuple | None:
        """Run the check on a specific position.

        Parameters
        ----------
        bases : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        tuple | None
            None if all nucleotide edits are sufficient, a tuple with
            error message otherwise.
        """
        for base in self._bases:
            if base != bases.reference and \
                    0 < bases[base] < self.min_edits_per_nucleotide:
                return (
                    'DISCARDING COLUMN edits={} < {}',
                    bases[base],
                    self.min_edits_per_nucleotide,
                )
        return None
