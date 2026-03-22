class RTResult(object):
    """RNA editing analysis for a single base position."""

    def __init__(self, bases, strand, contig, position):
        """
        RNA editing analysis for a single base position.

        Parameters:
            bases (compiledPosition): Bases found by REDItools
            strand (str): Strand of the position
            contig (str): Contig name
            position (int): Genomic position
        """
        self.contig = contig
        self.position = position + 1
        self.bases = bases
        self.strand = strand
        self._variants = bases.get_variants()

    @property
    def variants(self):
        """
        The detected variants at this position.

        Returns:
            list
        """
        ref = self.bases.ref
        return [f'{ref}{base}' for base in self._variants]

    @property
    def mean_quality(self):
        """
        Mean read quality of the base position.

        Returns:
            int
        """
        if self.bases:
            return sum(self.bases.qualities) / len(self.bases)
        return 0

    @property
    def edit_ratio(self):
        """
        Edit ratio.

        Returns:
            float
        """
        if self._variants:
            max_edits = max(self.bases[base] for base in self._variants)
            return max_edits / (max_edits + self.bases['REF'])
        return 0

    @property
    def reference(self):
        """
        Base in the reference genome.

        Returns:
            str
        """
        return self.bases.ref

    @property
    def depth(self):
        """
        How many reads cover the position. (post filtering).

        Returns:
            int
        """
        return len(self.bases)

    @property
    def per_base_depth(self):
        """
        How many reads had each base for this position.

        Returns:
            list
        """
        return list(iter(self.bases))
