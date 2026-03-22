from .qccheck import QCCheck

class CheckExclusions(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check if the bases object is in an excluded position.

        Parameters:
            rtools (REDItools): Object running the analysis
            bases (CompiledPosition): Data for analysis

        Returns:
            (bool): True if the position is not excluded
        """
        in_exclusions = rtools.exclude_regions.contains(
            bases.contig,
            bases.position,
        )
        if in_exclusions:
            self._log(rtools, 'DISCARD COLUMN in excluded region')
            return False
        return True
