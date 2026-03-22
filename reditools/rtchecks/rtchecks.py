"""Quality control for REDItools analyses."""

from .check_is_none import CheckIsNone

class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self):
        """Create a RTChecks object."""
        self.check_list = [CheckIsNone]

    def add(self, qc_check_class):
        """
        Add a QC check.

        Parameters:
            qc_check_class (QCCheck Sub Class): The check to perform
        """
        self.check_list.append(qc_check_class)

    def discard(self, qc_check_class):
        """
        Remove a QC check.

        Parameters:
            qc_check_class (QCCheck Sub Class): The check to discard
        """
        if qc_check_class in self.check_list:
            self.check_list.remove(qc_check_class)

    def check(self, rtools, bases):
        """
        Perform QC.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True of all checks pass, else false
        """
        for qc_check_class in self.check_list:
            if not qc_check_class.run_check(rtools, bases):
                return False
        return True
