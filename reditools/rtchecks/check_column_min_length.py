from .qccheck import QCCheck

class CheckColumnMinLength(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check read depth.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if the read depth is sufficient
        """
        if len(bases) < rtools.min_column_length:
            self._log(
                rtools,
                'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
                len(bases),
                rtools.min_column_length,
            )
            return False
        return True
