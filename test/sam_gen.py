"""Classes for SAM and FASTA file generation."""
from __future__ import annotations

import random
import re
from dataclasses import InitVar, dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from test.aligner import Aligner
from typing import Any, Iterator

from pysam import samtools


class Genome:
    """Genomic sequences object."""

    def __init__(self) -> None:
        """Initialize self."""
        self.contigs: dict[str, str] = {}

    def __getitem__(self, contig_name: str) -> str:
        """Retrive chromsomal sequence.

        Parameters
        ----------
        contig_name : str
            Chromosome name.

        Returns
        -------
        str
            Nucleotide sequence.
        """
        return self.contigs[contig_name]

    def add_contig(
        self,
        name: str | None=None,
        length: int=120,
        sequence: str | None=None,
    ) -> str:
        """Create a new chromosome.

        Parameters
        ----------
        name : str | None
            Chromosome name. If None, a name will be generated.
        length : int
            Sequence length (ignored if sequence is not None).
        sequence : str | None
            Nucleotide sequence. If None, a random sequence will be generated.

        Returns
        -------
        str
            Chromosome name.
        """
        if name is None:
            n_contigs = len(self.contigs)
            name = f"contig{len(self.contigs)}"
            while name in self.contigs:
                n_contigs += 1
                name = f"contig{len(self.contigs)}"
        if sequence is None:
            self.contigs[name] = self._random_seq(length)
        else:
            self.contigs[name] = sequence
        return name

    def save_to_fasta(self, filename: str) -> None:
        """Save the genome to a FASTA file.

        Parameters
        ----------
        filename : str
            Path to save file to.
        """
        with Path(filename).open("w") as stream:
            stream.writelines((
                f">{name} {idx}\n{sequence}\n"
                for idx, (name, sequence) in enumerate(self.contigs.items(), 1)
            ))
        samtools.faidx(filename)

    @classmethod
    def _random_seq(cls, length: int) -> str:
        """Generate a random nucleotide sequence.

        Parameters
        ----------
        length : int
            Sequence length.

        Returns
        -------
        str
            Random nucleotide sequence.
        """
        return "".join([random.choice("ACTG") for _ in range(length)])


@dataclass
class Sequence:
    """SAM entry.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.
    start : int
        Genomic start.
    flag : int
        SAM flags
    phred : list[int] | None
        PHRED scores (defaults to 30)
    mapq : int
        MAPQ score
    qname : str
        Read name.
    pnext : int
        Paired read start.
    """

    seq: str
    start: int
    flag: int = 0
    phred_list: InitVar[list | None] = None
    mapq: int = 255
    _cigar_str: str | None = None
    read_name: InitVar[str | None] = None
    pnext: int = 0

    read_n = 0
    flag_reverse_strand = 16
    phred_default = 30

    def __post_init__(
        self,
        phred_list: list[int] | None,
        read_name: str | None,
    ) -> None:
        """Post initialization.

        Parameters
        ----------
        phred : list[int] | None
            PHRED scores (defaults to 30)
        qname :  str | None
            Read name (one will be generated if None)
        """
        if phred_list is None:
            self.phred = [self.phred_default for _ in range(len(self.seq))]
        else:
            self.phred = phred_list

        if read_name is None:
            self.qname = self.next_read_name()
        else:
            self.qname = read_name

    def __len__(self) -> int:
        """Return sequence length.

        Returns
        -------
        int
            Sequence length.
        """
        return len(self.seq)

    def __str__(self) -> str:
        """Retrive nucleotide sequence.

        Returns
        -------
        str
            Nucleotide sequence.
        """
        return self.seq

    def tlen(self, ref_seq: str) -> int:
        """Calculate the transcript length.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.

        Returns
        -------
        int
            Transcript length.
        """
        cigar = self.cigar_str(ref_seq)
        tlen = 0
        for count, op in re.findall(r"(?P<count>\d+)(?P<op>[A-Z])", cigar):
            if op not in ("S", "I"):
                tlen += int(count)
        if self.flag & Sequence.flag_reverse_strand:
            return -tlen
        return tlen

    def cigar_str(self, ref_seq: str) -> str:
        """Generate alignment CIGAR string.

        Parameters
        ----------
        ref_seq : str
            Reference sequence

        Returns
        -------
        str
            CIGAR string.
        """
        if self._cigar_str is not None:
            return self._cigar_str
        alignment = Aligner().align(
            ref_seq[self.start:self.start + len(self)],
            str(self),
        )
        cigar_iter = self.assemble_cigar_list(*alignment)
        cigar_pieces = [f"{length}{op}" for length, op in cigar_iter]
        self._cigar_str = "".join(cigar_pieces)
        return self._cigar_str

    def make_pair(self) -> Sequence:
        """Generate SAM paired entry.

        Returns
        -------
        Sequence
            SAM paired sequence.
        """
        return Sequence(
            seq=self.seq,
            start=self.start,
            flag=self.pair_flag(self.flag),
            phred_list=self.phred,
            mapq=self.mapq,
            _cigar_str=self._cigar_str,
            read_name=self.qname,
            pnext=self.start,
        )

    @classmethod
    def pair_flag(cls, flag_value: int) -> int:
        """Create flags for paired read.

        Parameters
        ----------
        flag_value : int
            Flags of read to pair

        Returns
        -------
        int
            Flags for read mate.
        """
        return flag_value ^ 240

    @classmethod
    def cigar_op(cls, ref_base: str, query_base: str) -> str:
        """Determine CIGAR op from alignment base.

        Parameters
        ----------
        ref_base : str
            Reference sequence alignment base.
        query_base : str
            Query sequence alignment base.

        Returns
        -------
        str
            CIGAR operator.
        """
        if ref_base == "-":
            return "I"
        if query_base == "-":
            return "D"
        if query_base == ref_base:
            return "M"
        return "X"

    @classmethod
    def assemble_cigar_list(
        cls,
        algn_ref: str,
        algn_query: str,
    ) -> Iterator[tuple[int, str]]:
        """Iterate over CIGAR operators.

        Parameters
        ----------
        algn_ref : str
            Reference sequence alignment string.
        algn_query : str
            Query sequence alignment string.

        Yields
        ------
        tuple[int, str]
            CIGAR operator length and string.
        """
        last_op = None
        op_n = 0
        for ref_base, query_base in zip(algn_ref, algn_query):
            cigar_op = cls.cigar_op(ref_base, query_base)
            if cigar_op == last_op:
                op_n += 1
            else:
                if last_op is not None:
                    yield (op_n, cigar_op)
                last_op = cigar_op
                op_n = 1
        yield (op_n, cigar_op)

    @classmethod
    def next_read_name(cls) -> str:
        """Get the next available read name.

        Returns
        -------
        str
            Next read name in the form "read#"
        """
        cls.read_n += 1
        return f"read{cls.read_n}"


class SAM:
    """SAM file object."""

    def __init__(self) -> None:
        """Initialize self."""
        self.genome = Genome()
        self.reads: dict[str, list[Sequence]] = {}

    def __getitem__(self, contig_name: str) -> list[Sequence]:
        """Retrive reads for a specific contig.

        Parameters
        ----------
        contig_name : str
            Chromosome name.

        Returns
        -------
        list[Sequence]
            Reads aligned to the chromosome.
        """
        return self.reads[contig_name]

    def header(self) -> str:
        """Generate SAM header.

        Returns
        -------
        str
            Header.
        """
        header = ["@HD\tVN:1.5"]
        for contig_name, seq in self.genome.contigs.items():
            header.append(f"@SQ\tSN:{contig_name}\tLN:{len(seq)}")
        header.append(
            "@RG\tID:1\tSM:1_AAAAA\tLB:default\tPU:xxx.1\tPL:ILLUMINA",
        )
        header.append("@PG\tID:reditools\tPN:reditools\tCL:gen_sam.py")
        return "\n".join(header)

    def add_contig(
        self,
        contig_name: str | None=None,
        length: int=120,
        sequence: str | None=None,
    ) -> str:
        """Add a chromosome to the genomic reference.

        Parameters
        ----------
        contig_name : str | None
            Chromosome name (generated if not provided).
        length : int
            Chromosome size (used if sequence is None).
        sequence : str | None
            Genomic sequence (random if not provided).

        Returns
        -------
        str
            Chromosome name.
        """
        contig_name = self.genome.add_contig(contig_name, length, sequence)
        self.reads[contig_name] = []
        return contig_name

    def add_read(self, contig_name: str, sequence_obj: Sequence) -> None:
        """Add a read.

        Parameters
        ----------
        contig_name : str
            Chromosome the read aligns to.
        sequence_obj : Sequence
            Read.
        """
        self.reads[contig_name].append(sequence_obj)

    def add_read_pair(self, contig_name: str, sequence_obj: Sequence) -> None:
        """Add a read and generate its mate.

        Parameters
        ----------
        contig_name : str
            Chromosome the read aligns to.
        sequence_obj : Sequence
            Read.
        """
        self.add_read(contig_name, sequence_obj)
        self.add_read(contig_name, sequence_obj.make_pair())

    def sam_entries(self) -> Iterator[str]:
        """Iterate over SAM entries.

        Yields
        ------
        SAM file lines.
        """
        for contig, reads in self.reads.items():
            ref_seq = self.genome[contig]
            for sequence in reads:
                yield "\t".join([str(_) for _ in (
                    sequence.qname,
                    sequence.flag,
                    contig,
                    sequence.start + 1,
                    sequence.mapq,
                    sequence.cigar_str(ref_seq),
                    ("*", "=")[sequence.flag & 1],
                    sequence.pnext + 1,
                    sequence.tlen(ref_seq),
                    str(sequence),
                    "".join([self._phred(_) for _ in sequence.phred]),
                )])

    def save_to_sam(self, bam_filename: str, genome_filename: str) -> None:
        """Save SAM data to file.

        Parameters
        ----------
        bam_filename : str
            Path to save to.
        genome_filename : str
            Path to reference FASTA file.
        """
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                dir=".",
                suffix=".sam",
        ) as stream:
            sam_filename = stream.name
            stream.write(self.header())
            stream.write("\n")
            stream.write("\n".join(self.sam_entries()))
        md_sam = samtools.calmd(
            sam_filename,
            genome_filename,
            catch_stdout=True,
        )
        with Path(sam_filename).open("w") as stream:
            stream.writelines(md_sam)
        samtools.sort("-o", bam_filename, sam_filename)
        samtools.index(bam_filename)
        Path(sam_filename).unlink()

    @classmethod
    def _phred(cls, int_value: int) -> str:
        """Convert PHRED score to character.

        Parameters
        ----------
        int_value : int
            PHRED score.

        Returns
        -------
        str
            PHRED character.
        """
        return chr(33 + int_value)

def ntf(*args: Any, **kwargs: Any) -> str:  # noqa: ANN401
    """Create a new temporary file.

    Parameters
    ----------
    *args : Any
        Positional arguments to NamedTemporaryFile.
    **kwargs : Any
        Named arguments to NamedTemporaryFile.

    Returns
    -------
    str
        Path to temporary file.
    """
    with NamedTemporaryFile(  # type: ignore[call-overload]
            *args,
            delete=False,
            mode="w",
            **kwargs,
    ) as stream:
        return stream.name
