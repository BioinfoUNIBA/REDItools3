from .qccheck import QCCheck

class CheckColumnEditFrequency(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check the number of edits at the site.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if there are sufficient edits.
        """
        edits_no = len(bases) - bases['REF']
        if edits_no < rtools.min_edits:
            self._log(
                rtools,
                'DISCARDING COLUMN edits={} < {}',
                edits_no,
                rtools.min_edits,
            )
            return False
        return True
