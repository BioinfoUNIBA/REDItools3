import argparse

from reditools.compiled_position import RTResult


class CheckMinReadDepth:
    def __init__(self, options: argparse.Namespace):
        self.min_read_depth = options.min_read_depth

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return options.min_read_depth > 1

    def run_check(self, bases: RTResult) -> None | tuple:
        if len(bases) < self.min_read_depth:
            return (
                'DISCARDING COLUMN {} [MIN_READ_DEPTH={}]',
                len(bases),
                self.min_read_depth,
            )
        return None
