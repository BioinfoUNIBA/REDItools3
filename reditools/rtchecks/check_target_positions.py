from .qccheck import QCCheck

class CheckTargetPositions(QCCheck):
    @staticmethod
    def run_check(rtools, bases):
        """
        Check if the bases object is in a target region.

        Parameters:
            rtools (REDItools): Object running the analysis
            bases (CompiledPosition): Data for analysis

        Returns:
            (bool): True if the position is in a target region
        """
        in_targets = rtools.target_regions.contains(
            bases.contig,
            bases.position,
        )
        if not in_targets:
            self._log(
                rtools,
                'DISCARD COLUMN not in target regions',
            )
        return in_targets
