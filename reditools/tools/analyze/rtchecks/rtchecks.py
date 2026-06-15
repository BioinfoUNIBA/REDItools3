import argparse

from reditools.compiled_position import RTResult
from reditools.tools.analyze import rtchecks


class RTChecks(object):
    """
    Manage and execute a suite of checks on RNA editing results.

    Parameters
    ----------
    options : argparse.Namespace
        Command-line options that determine which checks are enabled.
    """

    def __init__(self, options: argparse.Namespace):
        """
        Initialize RTChecks with enabled check instances.

        Parameters
        ----------
        options : argparse.Namespace
            Command-line options used to filter and configure checks.
        """
        self.check_list = []

        for check in (
                rtchecks.CheckColumnEditFrequency,
                rtchecks.CheckColumnMinEdits,
                rtchecks.CheckMinReadDepth,
                rtchecks.CheckExclusions,
                rtchecks.CheckMaxEditingNucleotides,
                rtchecks.CheckTargetPositions,
                rtchecks.CheckVariants,
        ):
            if check.is_needed(options):
                self.check_list.append(check(options))

    def check(self, rtresult: RTResult) -> None | tuple:
        """
        Run all enabled checks against a set of base results.

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
