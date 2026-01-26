"""Quality control for REDItools analyses."""

from reditools import utils
from reditools.logger import Logger


class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self):
        """Create a RTChecks object."""
        self.check_list = [self.check_is_none]

    def add(self, function):
        """
        Add a QC check.

        Parameters:
            function (RTChecks method): The check to perform
        """
        self.check_list.append(function)

    def discard(self, function):
        """
        Remove a QC check.

        Parameters:
            function (RTChecks method): The check to discard
        """
        if function in self.check_list:
            self.check_list.remove(function)

    def check(self, rtools, bases):
        """
        Perform QC.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True of all checks pass, else false
        """
        return utils.check_list(
            self.check_list,
            bases=bases,
            rtools=rtools,
        )

    def check_column_min_length(self, rtools, bases):
        """
        Check read depth.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if the read depth is sufficient
        """
        if len(bases) < rtools.min_column_length:
            rtools.log(
                Logger.debug_level,
                'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
                len(bases),
                rtools.min_column_length,
            )
            return False
        return True

    # Really shouldn't use this one. I have to compute mean_q anyway
    def check_column_quality(self, rtools, bases):
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
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN mean_quality={} < {}',
                mean_q,
                rtools.min_read_quality,
            )
            return False
        return True

    def check_column_edit_frequency(self, rtools, bases):
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
            rtools.log(
                Logger.debug_level,
                'DISCARDING COLUMN edits={} < {}',
                edits_no,
                rtools.min_edits,
            )
            return False
        return True

    def check_column_min_edits(self, rtools, bases):
        """
        Check that there are sufficient edit events for each base.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if there are sufficient edits
        """
        for base in "ATCG":
            if base == bases.ref:
                continue
            if bases[base] < rtools.min_edits_per_nucleotide:
                rtools.log(
                    Logger.debug_level,
                    'DISCARDING COLUMN edits={} < {}',
                    bases[base],
                    rtools.min_edits_per_nucleotide,
                )
                return False
        return True

    def check_multiple_alts(self, bases, rtools):
        """
        Check that there is, at most, one alternate base.

        Parameters:
            bases (CompiledPosition): Base position under analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if there is zero or one alt
        """
        alts = bases.get_variants()
        if len(alts) < 2:
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN alts={} > 1',
                len(alts),
            )
            return False
        return True

    def check_is_none(self, bases, rtools):
        """
        Check if the bases object is None.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if bases is not None
        """
        if bases is None:
            rtools.log(Logger.debug_level, 'DISCARD COLUMN no reads')
            return False
        return True

    def check_target_positions(self, bases, rtools):
        """
        Check if the bases object is in a target region.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if the position is in a target region
        """
        in_targets = rtools.target_regions.contains(
            bases.contig,
            bases.position,
        )
        if not in_targets:
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN not in target regions',
            )
        return in_targets

    def check_exclusions(self, bases, rtools):
        """
        Check if the bases object is in an excluded position.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if the position is not excluded
        """
        in_exclusions = rtools.exclude_regions.contains(
            bases.contig,
            bases.position,
        )
        if in_exclusions:
            rtools.log(Logger.debug_level, 'DISCARD COLUMN in excluded region')
            return False
        return True

    def check_max_alts(self, bases, rtools):
        """
        Check that there are no more than a max number of alts.

        Parameters:
            bases (CompiledPosition): Base position under analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if there are n or fewer alts
        """

        alts = bases.get_variants()
        if len(alts) > rtools.max_alts:
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN alts={} > {}',
                len(alts),
                rtools.max_alts,
            )
            return False
        return True
