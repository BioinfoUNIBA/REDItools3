"""Organizational structure for tracking base coverage of genomic positions."""

from reditools.compiled_position import CompiledPosition


class CompiledReads:
    """Manager for CompiledPositions."""

    _strands = ('-', '+', '*')

    def __init__(
        self,
        strand=0,
        min_base_position=0,
        max_base_position=0,
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
            self.get_strand = lambda _: 2
        else:
            if strand == 1:
                self.forward_flags = {0, 99, 147}
            else:
                self.forward_flags = {16, 83, 163}
            self.get_strand = lambda _: _.flag in self.forward_flags

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
        return (_[2].upper() for _ in read.get_aligned_pairs(
            with_seq=True,
            matches_only=True,
        ))

    def _get_ref_from_fasta(self, read):
        pairs = read.get_aligned_pairs(matches_only=True)
        indices = [ref for _, ref in pairs]
        return self._ref.get_base(read.reference_name, *indices)

    def _prep_read(self, read):
        pairs = read.get_aligned_pairs(matches_only=True)
        for (read_pos, ref_pos), ref_base in zip(pairs, self._ref_seq(read)):
            if read_pos > read.query_length - self._qc['max_base_position']:
                break
            if read_pos < self._qc['min_base_position']:
                continue
            read_base = read.query_sequence[read_pos]
            if ref_base == 'N' or read_base == 'N':
                continue
            phred = read.query_qualities[read_pos]
            if phred < self._qc['min_base_quality']:
                continue
            yield (ref_pos, read_base, phred, ref_base)
