import argparse

from reditools.compiled_position import RTResult


class CheckColumnEditFrequency:
    def __init__(self, options: argparse.Namespace):
        self.min_edits = options.min_edits
    
    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return options.min_edits > 0

    def run_check(self, bases: RTResult) -> None | tuple:
        edits_no = len(bases) - bases['REF']
        if edits_no < self.min_edits:
            return (
                'DISCARDING COLUMN edits={} < {}',
                edits_no,
                self.min_edits,
            )
        return None
