from reditools import file_utils
from reditools.region_collection import RegionCollection


class CheckTargetPositions:
    def __init__(self, options):
        self.regions = RegionCollection()
        self.regions.add_regions(file_utils.read_bed_file(*options.bed_file))

    @classmethod
    def is_needed(cls, options):
        return options.bed_file is not None

    def run_check(self, bases):
        if not self.regions.contains(bases.contig, bases.position):
            return ('DISCARD COLUMN not in target regions',)
