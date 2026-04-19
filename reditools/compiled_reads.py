from typing import Iterator, Optional

from pysam import AlignedSegment

from reditools.compiled_position import CompiledPosition
from reditools.fasta_file import RTFastaFile


class RefFetch:
    """Helper class to fetch reference sequences.

    Depending on whether a FASTA file path is provided, it can either
    fetch the reference sequence from the FASTA file or directly from
    the AlignedSegment if MD tags are available.
    """

    def __init__(self, fasta_file_path: Optional[str] = None):
        """Initialize RefFetch.

        Parameters
        ----------
        fasta_file_path : str, optional
            Path to the reference FASTA file. If None, fetch from read.
        """
        if fasta_file_path:
            self.fasta_file = RTFastaFile(fasta_file_path)
            self._refseq_fn = self.get_ref_from_fasta
        else:
            self._refseq_fn = self.get_ref_from_read

    def get_refseq(self, read: AlignedSegment) -> Iterator[str]:
        """Fetch reference sequence. If a FASTA file was provided in the
        constructor, this function calla get_ref_from_fasta. Otherwise it
        calls get_ref_from_read.

        Parameters
        ----------
        read : AlignedSegment
            The alignment segment to extract the reference from.

        Returns
        -------
        Iterator[str]
            An iterator of reference nucleotides.
        """
        return self._refseq_fn(read)


    def get_ref_from_read(self, read: AlignedSegment) -> Iterator[str]:
        """Fetch reference sequence from the read itself.

        Parameters
        ----------
        read : AlignedSegment
            The alignment segment to extract the reference from.

        Returns
        -------
        Iterator[str]
            An iterator of reference nucleotides.
        """
        return (_[2].upper() for _ in read.get_aligned_pairs(
            with_seq=True,
            matches_only=True,
        ))

    def get_ref_from_fasta(self, read: AlignedSegment) -> Iterator[str]:
        """Fetch reference sequence from a FASTA file.

        Parameters
        ----------
        read : AlignedSegment
            The alignment segment providing the coordinates.

        Returns
        -------
        Iterator[str]
            An iterator of reference nucleotides.
        """
        pairs = read.get_aligned_pairs(matches_only=True)
        indices = [ref for _, ref in pairs]
        return self.fasta_file.get_base(
            read.reference_name,  # type: ignore
            *indices,
        )


class CompiledReads:
    """Class to compile read information into positions.

    It processes AlignedSegments and organizes base information by position
    and strand.
    """

    _strands = ('-', '+', '*')

    def __init__(
        self,
        strand: int = 0,
        min_base_position: int = 0,
        max_base_position: int = 0,
        min_base_quality: int = 0,
        fasta_file: Optional[str] = None,
    ):
        """Initialize CompiledReads.

        Parameters
        ----------
        strand : int, default 0
            Strandness of the data. 0 for unstranded, 1 for forward, 2 for
            reverse.
        min_base_position : int, default 0
            Minimum position from the start of the read to consider a base.
        max_base_position : int, default 0
            Minimum position from the end of the read to consider a base.
        min_base_quality : int, default 0
            Minimum Phred quality score to consider a base.
        fasta_file : str, optional
            Path to the reference FASTA file.
        """
        self._nucleotides: dict[int, CompiledPosition] = {}
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

    def add_reads(self, reads: list[AlignedSegment]) -> None:
        """Process and add reads to the compilation.

        Parameters
        ----------
        reads : list[AlignedSegment]
            List of AlignedSegment objects to process.
        """
        for read in reads:
            strand = self._strands[self.get_strand(read)]
            for pos, base, quality, ref in self._prep_read(read):
                if pos not in self._nucleotides:
                    self._nucleotides[pos] = CompiledPosition(
                        ref=ref,
                        position=pos,
                        contig=read.reference_name,  # type: ignore
                    )
                self._nucleotides[pos].add_base(quality, strand, base)

    def pop_range(self, start, stop) -> Iterator[CompiledPosition]:
        """Yield and remove CompiledPosition objects within a range.

        Parameters
        ----------
        start : int
            Start position (inclusive).
        stop : int
            Stop position (exclusive).

        Yields
        ------
        CompiledPosition
            The compiled position object.
        """
        for position in range(start, stop):
            if not self._nucleotides:
                break
            bases = self._nucleotides.pop(position, None)
            if bases is not None:
                yield bases

    def _prep_read(  # noqa: WPS231
            self,
            read: AlignedSegment,
    ) -> Iterator[tuple[int, str, int, str]]:
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
            read_base = read.query_sequence[read_pos]  # type: ignore
            if ref_base == 'N' or read_base == 'N':
                continue
            phred = read.query_qualities[read_pos]  # type: ignore
            if phred < self._qc['min_base_quality']:
                continue
            yield (ref_pos, read_base, phred, ref_base)

    def _unstranded_strand(self, read: AlignedSegment) -> int:
        return 2

    def _stranded_strand(self, read: AlignedSegment) -> int:
        return read.flag in self.forward_flags
