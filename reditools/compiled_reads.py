"""Organizational structure for tracking base coverage of genomic positions."""

from reditools.compiled_position import CompiledPosition
from reditools.fasta_file import RTFastaFile


class RefFetch:
    def __init__(self, fasta_file_path=None):
        if fasta_file_path:
            self.fasta_file = RTFastaFile(fasta_file_path)
            self.get_refseq = self.get_ref_from_fasta
        else:
            self.get_refseq = self.get_ref_from_read

    def get_ref_from_read(self, read):
        return (_[2].upper() for _ in read.get_aligned_pairs(
            with_seq=True,
            matches_only=True,
        ))

    def get_ref_from_fasta(self, read):
        pairs = read.get_aligned_pairs(matches_only=True)
        indices = [ref for _, ref in pairs]
        return self.fasta_file.get_base(read.reference_name, *indices)


class CompiledReads:
    """Manager for CompiledPositions."""

    _strands = ('-', '+', '*')

    def __init__(
        self,
        strand=0,
        min_base_position=0,
        max_base_position=0,
        min_base_quality=0,
        fasta_file=None,
    ):
        """
        Create a new CompiledReads object.

        Parameters:
            strand (int): Strand detection mode
            min_base_position (int): Left trims bases
            max_base_position (int): Right trims bases
            min_base_quality (int): Minimum base quality to report
            fasat_file (RTFastaFile): Optional genomic reference
        """
        self._nucleotides = {}
        if strand == 0:
            self.get_strand = self._unstranded_strand
        else:
            if strand == 1:
                self.forward_flags = {0, 99, 147}
            else:
                self.forward_flags = {16, 83, 163}
            self.get_strand = self._stranded_strand

        self.reference = RefFetch(fasta_file)

        self._qc = {
            'min_base_quality': min_base_quality,
            'min_base_position': min_base_position,
            'max_base_position': max_base_position,
        }

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
                if pos not in self._nucleotides:
                    self._nucleotides[pos] = CompiledPosition(
                        ref=ref,
                        position=pos,
                        contig=read.reference_name,
                    )
                self._nucleotides[pos].add_base(quality, strand, base)

    def pop_range(self, start, stop):
        """
        Iteratively calls pop across a genomic range.

        Params:
            start (int): Genomic start (zero index, inclusive)
            stop (int): Genomic stop (zero index, exclusive)

        Yields:
            CompiledPosition
        """
        for position in range(start, stop):
            if not self._nucleotides:
                break
            bases = self._nucleotides.pop(position, None)
            if bases is not None:
                yield bases

    def _prep_read(self, read):  # noqa: WPS231
        for (read_pos, ref_pos), ref_base in zip(
                read.get_aligned_pairs(matches_only=True),
                self.reference.get_refseq(read),
        ):
            # Right end trim
            if read_pos > read.query_length - self._qc['max_base_position']:
                break
            # Left end trim
            if read_pos < self._qc['min_base_position']:
                continue
            read_base = read.query_sequence[read_pos]
            if ref_base == 'N' or read_base == 'N':
                continue
            phred = read.query_qualities[read_pos]
            if phred < self._qc['min_base_quality']:
                continue
            yield (ref_pos, read_base, phred, ref_base)

    def _unstranded_strand(self, read):
        return 2

    def _stranded_strand(self, read):
        return read.flag in self.forward_flags
