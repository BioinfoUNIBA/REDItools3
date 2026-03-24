from dataclasses import dataclass, field
from reditools.compiled_position import CompiledPosition


@dataclass(slots=True)
class RTResult(object):
    """RNA editing analysis for a single base position."""
    bases: CompiledPosition
    strand: str
    contig: str
    position: int

    @property
    def variants(self):
        """
        The detected variants at this position.

        Returns:
            list
        """
        ref = self.bases.ref
        return [f'{ref}{base}' for base in self.bases.get_variants()]

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
        variants = self.variants
        if variants:
            max_edits = max(self.bases[base] for base in variants)
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
