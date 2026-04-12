
class CheckColumnEditFrequency:
    def __init__(self, options):
        self.min_edits = options.min_edits
    
    @classmethod
    def is_needed(cls, options):
        return options.min_edits > 0

    def run_check(self, bases):
        edits_no = len(bases) - bases['REF']
        if edits_no < self.min_edits:
            return (
                'DISCARDING COLUMN edits={} < {}',
                edits_no,
                self.min_edits,
            )
