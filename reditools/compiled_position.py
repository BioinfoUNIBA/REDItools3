"""Organizational structure for tracking base coverage of genomic positions."""
from dataclasses import dataclass, field


class RTResult:
    _base_order = 'ACGT'
    def __init__(self, compiled_position, strand):
        self.cp = compiled_position
        self.strand = strand

        self.reference = self.cp.ref
        self.position = self.cp.position
        self.contig = self.cp.contig

        self.counter = {_: 0 for _ in self._base_order}
        for base in self.cp.bases:
            self.counter[base] += 1

        self.variants = [
            f'{self.reference}{_}' for _ in self._base_order
            if self[_] and _ != self.reference
        ]


    def __getitem__(self, base):
        if base.upper() == 'REF':
            return self.counter[self.reference]
        return self.counter[base]

    def __iter__(self):
        return (self[base] for base in self._base_order)

    def __len__(self):
        return len(self.cp)

    def __getattr__(self, attr):
        return getattr(self.cp, attr)

    @property
    def edit_ratio(self):
        max_edits = 0
        for base, count in zip(self._base_order, self):
            if base != self.reference and count > max_edits:
                max_edits = count
        try:
            return max_edits / (self['REF'] + max_edits)
        except ZeroDivisionError:
            return 0

    @property
    def mean_quality(self):
        if len(self) == 0:
            return 0
        return sum(self.qualities) / len(self)

@dataclass
class CompiledPosition:
    ref: str
    position: int
    contig: str
    qualities: list[int] = field(default_factory=list) 
    strands: list[str] = field(default_factory=list) 
    bases: list[str] = field(default_factory=list) 

    _comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __len__(self):
        return len(self.bases)

    def add_base(self, quality, strand, base):
        self.bases.append(base)
        self.strands.append(strand)
        self.qualities.append(quality)

    def calculate_strand(self, threshold=0):
        pos_count = 0
        neg_count = 0
        for strand in self.strands:
            if strand == '+':
                pos_count += 1
            elif strand == '-':
                neg_count += 1
        if pos_count == neg_count:
            return '*'
        if pos_count / (pos_count + neg_count) >= threshold:
            return '+'
        if neg_count / (pos_count + neg_count) >= threshold:
            return '-'
        return '*'

    def filter_by_strand(self, strand):
        if strand == '*':
            return
        keep = range(len(self.bases))
        keep = [idx for idx in keep if self.strands[idx] == strand]
        self.qualities = [self.qualities[_] for _ in keep]
        self.strands = [self.strands[_] for _ in keep]
        self.bases = [self.bases[_] for _ in keep]

    def complement(self):
        self.bases = [self._comp[base] for base in self.bases]
        self.ref = self._comp[self.ref]
