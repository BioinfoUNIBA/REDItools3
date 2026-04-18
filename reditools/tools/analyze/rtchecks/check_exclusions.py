import argparse

from reditools import file_utils
from reditools.compiled_position import RTResult
from reditools.region_collection import RegionCollection


class CheckExclusions:
    def __init__(self, options: argparse.Namespace):
        self.regions = RegionCollection()
        self.regions.add_regions(file_utils.read_bed_file(*options.exclude_regions))

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return options.exclude_regions is not None

    def run_check(self, bases: RTResult) -> None | tuple:
        if self.regions.contains(bases.contig, bases.position):
            return ('DISCARD COLUMN in excluded region',)
        return None
