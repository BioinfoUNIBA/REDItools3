
class CheckMaxEditingNucleotides:
    def __init__(self, options):
        self.max_editing_nucleotides = options.max_editing_nucleotides

    @classmethod
    def is_needed(cls, options):
        return options.max_editing_nucleotides < 3

    def run_check(self, bases):
        variants = bases.variants
        if len(variants) > self.max_editing_nucleotides:
            return (
                'DISCARD COLUMN variants={} > {}',
                len(variants),
                self.max_editing_nucleotides,
            )

