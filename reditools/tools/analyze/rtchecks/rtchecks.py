"""Quality control for REDItools analyses."""
from reditools.tools.analyze import rtchecks


class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self, options):
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

    def check(self, bases):
        """
        Perform QC.

        Parameters:
            bases (CompiledPosition): Base position under analysis

        Returns:
            None if QC passed, else debug message (tuple)
        """
        for check_class in self.check_list:
            check_result = check_class.run_check(bases)
            if check_result is not None:
                return check_result
