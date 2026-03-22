from .qccheck import QCCheck

class CheckMaxAlts(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check that there are no more than a max number of alts.

        Parameters:
            rtools (REDItools): Object running the analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if there are n or fewer alts
        """

        alts = bases.get_variants()
        if len(alts) > rtools.max_alts:
            self._log(
                rtools,
                'DISCARD COLUMN alts={} > {}',
                len(alts),
                rtools.max_alts,
            )
            return False
        return True
