"""Organizational structure for tracking base coverage of genomic positions."""


class CompiledPosition(object):
    """Tracks base frequency for a genomic position."""

    _bases = 'ACGT'
    _comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, ref, contig, position):
        """
        Create a new compiled position.

        Parameters:
            ref (string): The reference base for this position
            contig (string): Chromosome name
            position (int): Genomic coordinate
        """
        self.qualities = []
        self.strands = []
        self.bases = []
        self.counter = False
        self.ref = ref
        self.contig = contig
        self.position = position

    def __len__(self):
        """
        Position depth.

        Returns:
            int
        """
        return len(self.bases)

    def __getitem__(self, base):
        """
        Frequency of a given nucleotide at this position.

        Parameters:
            base (str): The nucleotide (A, C, G, T, or REF)

        Returns:
            int: The total number of reads with the given base
        """
        if not self.counter:
            self.counter = {base: 0 for base in self._bases}
            for base_member in self.bases:
                self.counter[base_member] += 1
        if base.upper() == 'REF':
            return self.counter[self.ref]
        return self.counter[base]

    def __iter__(self):
        """
        Iterate over each base frequency.

        Returns:
            iterator
        """
        return (self[base] for base in self._bases)

    def add_base(self, quality, strand, base):
        """
        Add details for a base at this position.

        Parameters:
            quality (int): The quality of the read
            strand (str): The strand the base is on (+, -, or *)
            base (str): The nucleotide at the position( A, C, G, or T)
        """
        self.qualities.append(quality)
        self.strands.append(strand)
        self.bases.append(base)
        self.counter = False

    def complement(self):
        """Modify all the summarized nucleotides to their complements."""
        self.bases = [self._comp[base] for base in self.bases]
        self.ref = self._comp[self.ref]
        if not self.counter:
            return
        complements = self._comp.items()
        self.counter = {sb: self.counter[bs] for bs, sb in complements}

    def get_variants(self):
        """
        List all detected variants.

        Returns:
            list
        """
        alts = set(self._bases) - {self.ref}
        return [base for base in alts if self[base]]

    def get_strand(self, threshold=0):
        """
        Determine the mean strandedness of a position.

        Parameters:
            threshold (int): Confidence minimum for strand identification

        Returns:
            '+', '-', or '*'
        """
        strand_counts = {'+': 0, '-': 0, '*': 0}
        for idx in self.strands:
            strand_counts[idx] += 1
        total = strand_counts['+'] + strand_counts['-']
        if total == 0:
            return '*'

        strand = max(strand_counts, key=strand_counts.get)
        if strand_counts[strand] / total >= threshold:
            return strand
        return '*'

    def filter_by_strand(self, strand):
        """
        Remove all bases not on the strand.

        Parameters:
            strand (str): Either +, -, or *
        """
        keep = range(len(self.bases))
        keep = [idx for idx in keep if self.strands[idx] == strand]
        self.qualities = self._filter(self.qualities, keep)
        self.strands = self._filter(self.strands, keep)
        self.bases = self._filter(self.bases, keep)
        self.counter = False

    def _filter(self, lst, indx):
        return [lst[idx] for idx in indx]
