"""Quality control for REDItools analyses."""
from reditools.region_collection import RegionCollection
from reditools import file_utils
from reditools.tools.analyze import rtchecks


class RTChecks(object):
    """Quality control for REDItools analyses."""

    def __init__(self, options):
        self.check_list = []

        variants = set()
        for alt in options.variants:
            alt = alt.upper()
            if alt == 'ALL':
                variants = set()
                break
            if len(alt) != 2 or alt[0] not in 'ATCG' or alt[1] not in 'ATCG':
                raise ValueError(f'Bad variant: {alt}')
            variants.add(alt)

        if variants:
            options.variants = variants
            self.check_list.append(rtchecks.check_variants)

        if options.bed_file:
            options.target_regions = RegionCollection()
            options.target_regions.add_regions(
                file_utils.read_bed_file(*options.bed_file),
            )
            self.check_list.append(rtchecks.check_target_positions)

        if options.exclude_regions:
            excluded_regions = RegionCollection()
            excluded_regions.add_regions(
                file_utils.read_bed_file(*options.exclude_regions),
            )
            options.exclude_regions = excluded_regions
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


    def check(self, bases):
        """
        Perform QC.

        Parameters:
            bases (CompiledPosition): Base position under analysis
            rtools (REDItools): for logging

        Returns:
            (bool): True of all checks pass, else false
        """
        for qc_check_fn in self.check_list:
            msg = qc_check_fn(self.namespace, bases)
            if msg is not None:
                return msg
