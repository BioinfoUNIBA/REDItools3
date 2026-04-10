"""Quality control for REDItools analyses."""
from reditools.region_collection import RegionCollection
from reditools import file_utils
from reditools.tools.analyze import rtchecks
import re

class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self, options):
        self.check_list = []

        options.variants = self.verify_alts(options.variants)

        if options.bed_file:
            options.target_regions = self.create_region_collection(
                options.bed_file,
            )
            self.check_list.append(rtchecks.check_target_positions)

        if options.exclude_regions:
            self.exclude_regions = self.create_region_collection(
                options.exclude_regions,
            )
            self.check_list.append(rtchecks.check_exclusions)

        if options.max_editing_nucleotides < 3:
            self.check_list.append(rtchecks.check_max_editing_nucleotides)

        if options.min_read_depth > 1:
            self.check_list.append(rtchecks.check_min_read_depth)

        if options.min_edits > 0:
            self.check_list.append(rtchecks.check_column_edit_frequency)

        if options.min_edits_per_nucleotide > 0:
            self.check_list.append(rtchecks.check_column_min_edits)

        if options.splicing_file:
            regions = file_utils.load_splicing_file(
                options.splicing_file,
                options.splicing_span,
            )
            if not options.exclude_regions:
                options.exclude_regions = RegionCollection()
            options.exclude_regions.add_regions(regions)

        self.namespace = options

    def verify_alts(self, variants):
        alt_pa = re.compile('[ACTG]{2}')
        for alt in variants:
            alt = alt.upper()
            if alt == 'ALL':
                return set()
            if not alt_pa.fullmatch(alt):
                raise ValueError(f'Bad variant: {alt}')
        return set(variants)

    def create_region_collection(self, file_list):
        rc = RegionCollection()
        rc.add_regions(file_utils.read_bed_file(*file_list))
        return rc

    def check(self, bases):
        """
        Perform QC.

        Parameters:
            bases (CompiledPosition): Base position under analysis

        Returns:
            None if QC passed, else debug message (tuple)
        """
        return next((_(self.namespace, bases) for _ in self.check_list), None)
