from reditools.region_collection import RegionCollection
from reditools import file_utils

class CheckExclusions:
    def __init__(self, options):
        self.regions = RegionCollection()
        self.regions.add_regions(file_utils.read_bed_file(*options.exclude_regions))

    @classmethod
    def is_needed(cls, options):
        return options.exclude_regions is not None

    def run_check(self, bases):
        if self.regions.contains(bases.contig, bases.position):
            return ('DISCARD COLUMN in excluded region',)
