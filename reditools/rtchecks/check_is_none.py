from .qccheck import QCCheck

class CheckIsNone(QCCheck):
    @staticmethod
    def run_Check(rtools, bases):
        """
        Check if the bases object is None.

        Parameters:
            rtools (REDItools): Object running the analysis
            bases (CompiledPosition): Data for analysis

        Returns:
            (bool): True if bases is not None
        """
        if bases is None:
            self._log(rtools, 'DISCARD COLUMN no reads')
            return False
        return True
