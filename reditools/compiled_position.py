from dataclasses import dataclass, field
from typing import Iterator


@dataclass
class CompiledPosition:
    """Class to store compiled information for a specific genomic position.

    Attributes
    ----------
    ref : str
        The reference nucleotide at this position.
    position : int
        The genomic position (0-based).
    contig : str
        The name of the contig/chromosome.
    qualities : list[int]
        List of Phred quality scores for each base at this position.
    strands : list[str]
        List of strand orientations ('+', '-', or '*') for each base.
    bases : list[str]
        List of nucleotides observed at this position.
    """

    ref: str
    position: int
    contig: str
    qualities: list[int] = field(default_factory=list)
    strands: list[str] = field(default_factory=list)
    bases: list[str] = field(default_factory=list)

    _comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __len__(self) -> int:
        """Return the number of bases at this position.

        Returns
        -------
        int
            Number of bases.
        """
        return len(self.bases)

    def add_base(self, quality: int, strand: str, base: str) -> None:
        """Add a base observation to this position.

        Parameters
        ----------
        quality : int
            Phred quality score of the base.
        strand : str
            Strand orientation of the read containing the base
            ('+', '-', or '*').
        base : str
            The observed nucleotide.
        """
        self.bases.append(base)
        self.strands.append(strand)
        self.qualities.append(quality)

    def calculate_strand(self, threshold: float = 0) -> str:
        """Calculate the consensus strand based on observations.

        Parameters
        ----------
        threshold : float, default 0
            The fraction of observations required to assign a strand.

        Returns
        -------
        str
            The calculated strand ('+', '-', or '*').
        """
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

    def filter_by_strand(self, strand: str) -> None:
        """Filter observations to keep only those from a specific strand.

        Parameters
        ----------
        strand : str
            The strand to keep ('+', '-', or '*'). If '*', no filtering is done.
        """
        if strand == '*':
            return
        keep = [
            idx for idx in range(len(self.bases))
            if self.strands[idx] == strand
        ]
        self.qualities = [self.qualities[_] for _ in keep]
        self.strands = [self.strands[_] for _ in keep]
        self.bases = [self.bases[_] for _ in keep]

    def complement(self) -> None:
        """Replace all bases and the reference with their complements."""
        self.bases = [self._comp[base] for base in self.bases]
        self.ref = self._comp[self.ref]


class RTResult:
    """Class to represent the results of REDItools analysis at a position.

    Attributes
    ----------
    cp : CompiledPosition
        The compiled position data.
    strand : str
        The strand assigned to this position.
    reference : str
        The reference nucleotide.
    position : int
        The genomic position.
    contig : str
        The contig name.
    counter : dict
        A dictionary counting occurrences of each nucleotide.
    variants : list
        A list of observed variants (e.g., ['AG']).
    """

    _base_order = 'ACGT'

    def __init__(self, compiled_position: CompiledPosition, strand: str):
        """Initialize RTResult.

        Parameters
        ----------
        compiled_position : CompiledPosition
            The compiled position data.
        strand : str
            The strand assigned to this position.
        """
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

    def __getitem__(self, base: str):
        """Get the count of a specific base or the reference base.

        Parameters
        ----------
        base : str
            The nucleotide ('A', 'C', 'G', 'T') or 'REF'.

        Returns
        -------
        int
            The count of the requested base.
        """
        if base.upper() == 'REF':
            return self.counter[self.reference]
        return self.counter[base]

    def __iter__(self) -> Iterator[int]:
        """Iterate over the counts of bases in ACGT order.

        Yields
        ------
        int
            The count of each base in order: 'A, 'C', 'G', 'T'.
        """
        return (self[base] for base in self._base_order)

    def __len__(self) -> int:
        """Return the total number of reads at this position.

        Returns
        -------
        int
            Total reads.
        """
        return len(self.cp)

    @property
    def edit_ratio(self) -> float:
        """Calculate the editing ratio at this position.

        The ratio of the most frequent variant to the sum of reference
        and that variant.

        Returns
        -------
        float
            The editing ratio.
        """
        max_edits = 0
        for base, count in zip(self._base_order, self):
            if base != self.reference and count > max_edits:
                max_edits = count
        try:
            return max_edits / (self['REF'] + max_edits)
        except ZeroDivisionError:
            return 0

    @property
    def mean_quality(self) -> float:
        """Calculate the mean quality score at this position.

        Returns
        -------
        float
            The mean quality score.
        """
        if len(self) == 0:
            return 0
        return sum(self.cp.qualities) / len(self)
