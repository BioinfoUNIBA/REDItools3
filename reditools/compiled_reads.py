"""Organizational structure for tracking base coverage of genomic positions."""

from reditools.compiled_position import CompiledPosition

inf = float('inf')


class CompiledReads(object):
    """Manager for CompiledPositions."""

    _strands = ('-', '+', '*')

    def __init__(
        self,
        strand=0,
        min_base_position=0,
        max_base_position=inf,
        min_base_quality=0,
    ):
        """
        Create a new CompiledReads object.

        Parameters:
            strand (int): Strand detection mode
            min_base_position (int): Left trims bases
            max_base_position (int): Right trims bases
            min_base_quality (int): Minimum base quality to report
        """
        self._nucleotides = {}
        if strand == 0:
            self.get_strand = lambda read: read.is_reverse
        else:
            self.get_strand = self._get_strand

        self._strand_one = strand == 1
        self._ref = None
        self._ref_seq = self._get_ref_from_read

        self._qc = {
            'min_base_quality': min_base_quality,
            'min_base_position': min_base_position,
            'max_base_position': max_base_position,
        }

    def add_reference(self, ref):
        """
        Add a reference FASTA file to use.

        Parameters:
            ref (RTFastaFile): Reference sequence
        """
        self._ref = ref
        self._ref_seq = self._get_ref_from_fasta

    def add_reads(self, reads):
        """
        Add iterable of pysam reads to the object.

        The reads are broken down. into individual nucleotides that are
        tracked by chromosomal location.

        Parameters:
            reads (iterable): pysam reads
        """
        for read in reads:
            strand = self._strands[self.get_strand(read)]
            for pos, base, quality, ref in self._prep_read(read):
                try:
                    self._nucleotides[pos].add_base(quality, strand, base)
                except KeyError:
                    self._nucleotides[pos] = CompiledPosition(
                        ref=ref,
                        position=pos,
                        contig=read.reference_name,
                    )
                    self._nucleotides[pos].add_base(quality, strand, base)

    def pop(self, position):
        """
        Remove and return the CompiledPosition at position.

        Method returns None if the position is empty.

        Parameters:
            position (int): The chromosomal location to pop

        Returns:
            A CompiledPosition or None if position is empty.
        """
        return self._nucleotides.pop(position, None)

    def is_empty(self):
        """
        Determine if there are any CompiledPositions still in the object.

        Returns:
            True if the object is empty, else False
        """
        return not self._nucleotides

    def _get_ref_from_read(self, read):
        return list(read.get_reference_sequence().upper())

    def _get_ref_from_fasta(self, read):
        pairs = read.get_aligned_pairs(matches_only=True)
        indices = [ref for _, ref in pairs]
        return self._ref.get_base(read.reference_name, *indices)

    def _qc_base_position(self, read, position):
        return read.query_length - position >= self._qc['max_base_position']

    def _prep_read(self, read):
        pairs = read.get_aligned_pairs(matches_only=True)
        seq = read.query_sequence.upper()
        qualities = read.query_qualities
        ref_seq = self._ref_seq(read)
        while pairs and pairs[0][0] < self._qc['min_base_position']:
            pairs.pop(0)
            ref_seq.pop(0)
        if not pairs:
            return

        while pairs and self._qc_base_position(read, pairs[0][0]):
            offset, ref_pos = pairs.pop(0)
            ref_base = ref_seq.pop(0)
            if ref_base != 'N' != seq[offset]:
                if qualities[offset] >= self._qc['min_base_quality']:
                    yield (ref_pos, seq[offset], qualities[offset], ref_base)

    def _get_strand(self, read):
        return read.is_read2 ^ self._strand_one ^ read.is_reverse
