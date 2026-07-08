"""Manage and execute a suites of checks and filters on RNA editing results."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools.tools.analyze import rtchecks

if TYPE_CHECKING:
    import argparse

    from reditools.compiled_position import RTResult

all_checks = (
    rtchecks.CheckColumnEditFrequency,
    rtchecks.CheckColumnMinEdits,
    rtchecks.CheckMinReadDepth,
    rtchecks.CheckExclusions,
    rtchecks.CheckMaxEditingNucleotides,
    rtchecks.CheckTargetPositions,
    rtchecks.CheckVariants,
)

class RTChecks:
    """Manage and execute a suite of checks on RNA editing results.

    Parameters
    ----------
    options : argparse.Namespace
        Command-line options that determine which checks are enabled.
    """

    def __init__(self, options: argparse.Namespace) -> None:
        """Initialize RTChecks with enabled check instances.

        Parameters
        ----------
        options : argparse.Namespace
            Command-line options used to filter and configure checks.
        """
        self.check_list = [
            check(options) for check in all_checks if check.is_needed(options)
        ]

    def check(self, rtresult: RTResult) -> None | tuple:
        """Run all enabled checks against a set of base results.

        Parameters
        ----------
        rtresult : RTResult
            The REDItools analysis result for a position.

        Returns
        -------
        Optional[tuple]
            The result of the first failing check, or None if all checks pass.
        """
        generator = (_.run_check(rtresult) for _ in self.check_list)
        return next((_ for _ in generator if _ is not None), None)
