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

    def check_splice_positions(self, rtools, bases):
        """
        Check if the contig and position are in a splice site.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis

        Returns:
            (bool): True if the position is not a splice site.
        """
        contig = bases.contig
        if bases.position in rtools.splice_positions.get(contig, []):
            rtools.log(
                Logger.debug_level,
                '[SPLICE_SITE] Discarding ({}, {}) because in splice site',
                contig,
                bases.position,
            )
            return False
        return True

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
        for num_edits in bases.get_min_edits():
            if 0 < num_edits < rtools.min_edits_per_nucleotide:
                rtools.log(
                    Logger.debug_level,
                    'DISCARDING COLUMN edits={} < {}',
                    num_edits,
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
        if bases.position not in rtools.target_positions.get(bases.contig, []):
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN not in target positions',
            )
            return False
        return True

    def check_ref(self, bases, rtools):
        """
        Check if the reference base is of interest.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if reference base was specified
        """
        if bases.ref not in rtools.include_refs:
            rtools.log(
                Logger.debug_level,
                'DISCARD COLUMN base "{}" not listed for reporting',
                bases.ref,
            )
            return False
        return True

    def check_exclusions(self, bases, rtools):
        """
        Check if the bases object is in an excluded position.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if the position is not excluded
        """
        if bases.position in rtools.exclude_positions.get(bases.contig, []):
            rtools.log(Logger.debug_level, 'DISCARD COLUMN in excluded region')
            return False
        return True

    def check_specific_edits(self, bases, rtools):
        """
        Check whether specified edits are present.

        Parameters:
            bases (CompiledPosition): Data for analysis
            rtools (REDItools): Object running the analysis

        Returns:
            (bool): True if the edit was specified
        """
        for ref, alt in rtools.specific_edits:
            if not bases[ref] or not bases[alt]:
                rtools.log(
                    Logger.debug_level,
                    'DISCARD COLUMN edit "{}" not specified for output',
                    ref + alt,
                )
                return False
        return True
