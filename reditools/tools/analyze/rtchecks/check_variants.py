import argparse
import re

from reditools.compiled_position import RTResult


class CheckVariants:
    def __init__(self, options: argparse.Namespace):
        pa = re.compile('[ATCG]{2}', re.IGNORECASE)
        bad_alt = next(
            (_ for _ in options.variants if not pa.fullmatch(_)),
            None,
        )
        if bad_alt is not None:
            raise ValueError(
                'Bad variant ({bad_alt}). Must be two bases (e.g. AG).'
            )
        self.variants = {_.upper() for _ in options.variants}

    @classmethod
    def is_needed(cls, options: argparse.Namespace) -> bool:
        return 'ALL' not in [_.upper() for _ in options.variants]

    def run_check(self, bases: RTResult) -> None | tuple:
        if any(_ in self.variants for _ in bases.variants):
            return None
        return (
            'DISCARD COLUMN Edits {} not in requested alts {}',
            bases.variants,
            self.variants,
        )
