
class CheckMaxEditingNucleotides:
    def __init__(self, options):
        self.max_editing_nucleotides = options.max_editing_nucleotides

    @classmethod
    def is_needed(cls, options):
        return options.max_editing_nucleotides < 3

    def run_check(self, bases):
        alts = bases.alts
        if len(alts) > self.max_editing_nucleotides:
            return (
                'DISCARD COLUMN alts={} > {}',
                len(alts),
                self.max_editing_nucleotides,
            )

