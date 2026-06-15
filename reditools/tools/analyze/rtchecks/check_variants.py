import argparse
import re

from reditools.compiled_position import RTResult


class CheckVariants:
    """
    Check if detected variants match specified allowed variants.

    Parameters
    ----------
    options : argparse.Namespace
        Command-line options containing allowed variants.
    """

    def __init__(self, options: argparse.Namespace):
        """
        Initialize CheckVariants with allowed variants.

        Parameters
        ----------
        options : argparse.Namespace
            Options containing 'variants', a list of two-base strings.

        Raises
        ------
        ValueError
            If a variant is not exactly two bases (e.g., 'AG').
        """
        pa = re.compile('[ATCG]{2}', re.IGNORECASE)
        bad_alt = next(
            (_ for _ in options.variants if not pa.fullmatch(_)),
            None,
        )
        if bad_alt is not None:
            raise ValueError(
                f'Bad variant ({bad_alt}). Must be two bases (e.g. AG).'
            )
        self.variants = {_.upper() for _ in options.variants}

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        """
        Determine if the variant check is required.

        Parameters
        ----------
        options : argparse.Namespace
            Options to check for variant requirements.

        Returns
        -------
        bool
            True if specific variants are required, False if 'ALL' is in the
            variant list.
        """
        return 'ALL' not in [_.upper() for _ in options.variants]

    def run_check(self, rtresult: RTResult) -> None | tuple:
        """
        Verify that detected variants are among the allowed ones.

        Parameters
        ----------
        rtresult : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        Optional[tuple]
            None if at least one variant is allowed, otherwise a failure
            message.
        """
        if any(_ in self.variants for _ in rtresult.variants):
            return None
        return (
            'DISCARD COLUMN Edits {} not in requested alts {}',
            rtresult.variants,
            self.variants,
        )
