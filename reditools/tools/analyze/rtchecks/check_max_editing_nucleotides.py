import argparse

from reditools.compiled_position import RTResult


class CheckMaxEditingNucleotides:
    """Check if a position has at most a certain number of editing nucleotides.

    Attributes
    ----------
    max_editing_nucleotides : int
        The maximum allowed number of editing nucleotides.
    """

    def __init__(self, options: argparse.Namespace):
        """Initialize CheckMaxEditingNucleotides.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options containing max_editing_nucleotides.
        """
        self.max_editing_nucleotides = options.max_editing_nucleotides

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
            True if max_editing_nucleotides < 3, False otherwise.
        """
        return options.max_editing_nucleotides < 3

    def run_check(self, rtresult: RTResult) -> None | tuple:
        """Run the check on a specific position.

        Parameters
        ----------
        rtresult : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        None | tuple
            None if number of variants is within limits, a tuple with
            error message otherwise.
        """
        variants = rtresult.variants
        if len(variants) > self.max_editing_nucleotides:
            return (
                'DISCARD COLUMN variants={} > {}',
                len(variants),
                self.max_editing_nucleotides,
            )
        return None

