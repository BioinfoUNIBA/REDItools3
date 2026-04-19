import argparse

from reditools.compiled_position import RTResult


class CheckMaxEditingNucleotides:
    def __init__(self, options: argparse.Namespace):
        self.max_editing_nucleotides = options.max_editing_nucleotides

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return options.max_editing_nucleotides < 3

    def run_check(self, bases: RTResult) -> None | tuple:
        variants = bases.variants
        if len(variants) > self.max_editing_nucleotides:
            return (
                'DISCARD COLUMN variants={} > {}',
                len(variants),
                self.max_editing_nucleotides,
            )
        return None

