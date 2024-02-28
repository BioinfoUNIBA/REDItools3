"""Quality control for REDItools analyses."""

from reditools import utils
from reditools.logger import Logger


class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self):
        """Create a RTChecks object."""
        self.check_list = set()

    def add(self, function):
        """
        Add a QC check.

        Parameters:
            function (RTChecks method): The check to perform
        """
        self.check_list.add(function)

    def discard(self, function):
        """
        Remove a QC check.

        Parameters:
            function (RTChecks method): The check to discard
        """
        self.check_list.discard(function)

    def check(self, rtools, bases, contig, position):
        """
        Perform QC.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True of all checks pass, else false
        """
        return utils.check_list(self.check_list, bases=bases)

    def check_splice_positions(self, rtools, bases, contig, position):
        """
        Check if the contig and position are in a splice site.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True if the position is not a splice site.
        """
        if position in rtools.splice_positions.get('contig', []):
            rtools.log(
                Logger.debug_level,
                '[SPLICE_SITE] Discarding ({}, {}) because in splice site',
                contig,
                position,
            )
            return False
        return True

    def check_poly_positions(self, rtools, bases, contig, position):
        """
        Check if the position is omoplymeric.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True if the position is not omopolymeric
        """
        if position in rtools.omopolymeric_positions.get(contig, []):
            rtools.log(
                Logger.debug_level,
                '[OMOPOLYMERIC] Discarding position ({}, {})' +
                'because omopolymeric',
                contig,
                position,
            )
            return False
        return True

    def check_column_min_length(self, rtools, bases, contig, position):
        """
        Check read depth.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True if the read depth is sufficient
        """
        if len(bases) < rtools.min_column_length:
            self._log(
                Logger.debug_level,
                'DISCARDING COLUMN {} [MIN_COLUMN_LEGNTH={}]',
                len(bases),
                rtools.min_column_length,
            )
            return False
        return True

    # Really shouldn't use this one. I have to compute mean_q anyway
    def check_column_quality(self, rtools, bases, contig, position):
        """
        Check mean quality of the position.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True if quality is sufficient
        """
        if bases:
            mean_q = sum(bases.qualities) / len(bases)
        else:
            mean_q = 0
        if mean_q < rtools.min_read_quality:
            self._log(
                Logger.debug_level,
                'DISCARD COLUMN mean_quality={} < {}',
                mean_q,
                rtools.min_read_quality,
            )
            return False
        return True

    def check_column_edit_frequency(self, rtools, bases, contig, position):
        """
        Check the number of edits at the site.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

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

    def check_column_min_edits(self, rtools, bases, contig, position):
        """
        Check that there are sufficient edit events for each base.

        Parameters:
            rtools (REDItools): Object performing analysis
            bases (CompiledPosition): Base position under analysis
            contig (str): Current contig
            position (int): Current position

        Returns:
            (bool): True if there are sufficient edits
        """
        for num_edits in bases.getmin_edits():
            if 0 < num_edits < rtools.min_edits_per_nucleotide:
                rtools.log(
                    Logger.debug_level,
                    'DISCARDING COLUMN edits={} < {}',
                    num_edits,
                    rtools.min_edits_per_nucleotide,
                )
                return False
        return True
