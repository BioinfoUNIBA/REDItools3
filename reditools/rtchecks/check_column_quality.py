from .qccheck import QCCheck

# Really shouldn't use this one. I have to compute mean_q anyway
class CheckColumnQuality(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check mean quality of the position.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if quality is sufficient
        """
        if bases:
            mean_q = sum(bases.qualities) / len(bases)
        else:
            mean_q = 0
        if mean_q < rtools.min_read_quality:
            self._log(
                rtools,
                'DISCARD COLUMN mean_quality={} < {}',
                mean_q,
                rtools.min_read_quality,
            )
            return False
        return True
