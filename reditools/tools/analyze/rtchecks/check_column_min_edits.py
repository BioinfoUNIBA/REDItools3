import argparse

from reditools.compiled_position import RTResult


class CheckColumnMinEdits:
    _bases = ('A', 'T', 'C', 'G')

    def __init__(self, options: argparse.Namespace):
        self.min_edits_per_nucleotide = options.min_edits_per_nucleotide
    
    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return options.min_edits_per_nucleotide > 0

    def run_check(self, bases: RTResult) -> tuple | None:
        for base in self._bases:
            if base != bases.reference and \
                    0 < bases[base] < self.min_edits_per_nucleotide:
                return (
                    'DISCARDING COLUMN edits={} < {}',
                    bases[base],
                    self.min_edits_per_nucleotide,
                )
        return None
