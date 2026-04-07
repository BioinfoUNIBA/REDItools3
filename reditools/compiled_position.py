"""Organizational structure for tracking base coverage of genomic positions."""


class CompiledPosition:
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
        self._strand = None

    def __len__(self):
        """
        Position depth.
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
        Iterate over each base frequency in order A, C, G, T.
        """
        return (self[base] for base in self._bases)

    def add_base(self, quality, strand, base):
        """
        Add details for a base at this position.

        Parameters:
            quality (int): The quality of the read
            strand (str): The strand the base is on (+, -, or *)
            base (str): The nucleotide at the position (A, C, G, or T)
        """
        self.qualities.append(quality)
        self.strands.append(strand)
        self.bases.append(base)
        self.counter = False

    def complement(self):
        """
        Modify all the summarized nucleotides to their complements.
        """
        self.bases = [self._comp[base] for base in self.bases]
        self.ref = self._comp[self.ref]
        if not self.counter:
            return
        complements = self._comp.items()
        self.counter = {sb: self.counter[bs] for bs, sb in complements}

    @property
    def reference(self):
        """
        str: Reference base at this position.
        (alias for ref)
        """
        return self.ref

    @property
    def alts(self):
        """
        list: Detected alternate bases.
        """
        return [b for b in self._bases if self[b] and b != self.ref]

    @property
    def variants(self):
        """
        list: Observed edits at this position (e.g. AG).
        """
        return [f'{self.ref}{base}' for base in self.alts]

    @property
    def mean_quality(self):
        """
        int: Mean read quality of the base position.
        """
        if len(self) == 0:
            return 0
        return sum(self.qualities) / len(self)

    @property
    def edit_ratio(self):
        """
        float: Edit ratio as most edited base frequency divided by sum of most
        edited base and reference base.
        """
        max_edits = 0
        for base, count in zip(self._bases, self):
            if base != self.ref and count > max_edits:
                max_edits = count
        try:
            return max_edits / (self['REF'] + max_edits)
        except ZeroDivisionError:
            return 0

    @property
    def depth(self):
        """
        int: Number of reads covering this position.
        (alias for __len__)
        """
        return len(self)

    @property
    def strand(self):
        """
        str: Get the strand information for this position ('+', '-', or '*')

        Raises:
            ValueError: If calculate_strand has not yet been run.
        """
        if self._strand is None:
            raise ValueError('Must run calculate_strand first.')
        return self._strand

    def calculate_strand(self, threshold=0):
        """
        Determine the strandedness.

        Parameters:
            threshold (int): Confidence minimum for strand identification

        Returns:
            '+', '-', or '*'
        """
        pos_count = 0
        neg_count = 0
        for strand in self.strands:
            if strand == '+':
                pos_count += 1
            elif strand == '-':
                neg_count += 1
        if pos_count == neg_count:
            self._strand = '*'
        elif pos_count / (pos_count + neg_count) >= threshold:
            self._strand = '+'
        elif neg_count / (pos_count + neg_count) >= threshold:
            self._strand = '-'
        else:
            self._strand = '*'
        return self._strand

    def filter_by_strand(self):
        """
        Remove all bases not on the strand.
        """
        if self.strand == '*':
            return
        keep = range(len(self.bases))
        keep = [idx for idx in keep if self.strands[idx] == self.strand]
        self.qualities = self._filter(self.qualities, keep)
        self.strands = self._filter(self.strands, keep)
        self.bases = self._filter(self.bases, keep)
        self.counter = False

    def _filter(self, lst, indx):
        return [lst[idx] for idx in indx]

    @property
    def per_base_depth(self):
        """
        list: Get the depth per base for this position in order A, C, G, T.
        """
        return list(self)
