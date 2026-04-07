"""Quality control for REDItools analyses."""


class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self):
        """Create a RTChecks object."""
        self.check_list = []

    def add(self, qc_check_fn):
        """
        Add a QC check.

        Parameters:
            qc_check_fn (function): The check to perform
        """
        self.check_list.append(qc_check_fn)

    def discard(self, qc_check_fn):
        """
        Remove a QC check.

        Parameters:
            qc_check_fn (function): The check to discard
        """
        if qc_check_fn in self.check_list:
            self.check_list.remove(qc_check_fn)

    def check(self, rtools, bases):
        """
        Perform QC.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True of all checks pass, else false
        """
        for qc_check_fn in self.check_list:
            if not qc_check_fn(rtools, bases):
                return False
        return True
