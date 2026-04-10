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
        return next((_.run_check(bases) for _ in self.check_list), None)
