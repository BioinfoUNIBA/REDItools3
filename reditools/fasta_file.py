from __future__ import annotations

from typing import TYPE_CHECKING

from pysam.libcfaidx import FastaFile as PysamFastaFile

if TYPE_CHECKING:
    from types import TracebackType
    from typing import Iterator


class MissingContigError(KeyError):
    def __init__(self, contig_name: str) -> None:
        self.message = f"Reference name {contig_name} not found in FASTA file."
        super().__init__(self.message)

class PastContigEndError(IndexError):
    def __init__(self, contig_name: str, position: int) -> None:
        self.message = (
            f"Base position {position} is outside the bounds of "
            "{contig}. Are you using the correct reference?"
        )
        super().__init__(self.message)

class RTFastaFile:
    """A wrapper around pysam.FastaFile for genomic sequence access."""

    def __init__(self, filename: str) -> None:
        """Initialize the RTFastaFile.

        Parameters
        ----------
        *args
            Arguments passed to pysam.FastaFile.
        **kwargs
            Keyword arguments passed to pysam.FastaFile.
        """
        self.pysam_fasta_file = PysamFastaFile(filename)

    def __enter__(self) -> RTFastaFile:
        return self

    def __exit__(
        self,
        typ: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Exit the runtime context related to this object.

        Parameters
        ----------
        typ : type[BaseException] | None
            The exception type.
        exc : BaseException | None
            The exception value.
        tb : TracebackType | None
            The traceback.
        """
        self.pysam_fasta_file.close()

    def get_base(self, contig: str, *position: int) -> Iterator[str]:
        """Retrieve bases at specified positions from a contig.

        Parameters
        ----------
        contig : str
            The name of the contig or chromosome.
        *position : int
            One or more 0-based positions to retrieve bases for.

        Returns
        -------
        Iterator[str]
            An iterator over the upper-case bases at the specified positions.

        Raises
        ------
        MissingContigError
            If the contig is not found in the FASTA file.
        PastContigEndError
            If a position is outside the bounds of the contig.
        """

        if contig not in self.pysam_fasta_file:
            if contig.startswith("chr"):
                new_contig = contig.replace("chr", "")
            else:
                new_contig = f"chr{contig}"
            if new_contig not in self.pysam_fasta_file:
                raise MissingContigError(contig)
            contig = new_contig
        sorted_pos = sorted(position)
        seq = self.pysam_fasta_file.fetch(
            contig,
            sorted_pos[0],
            sorted_pos[-1] + 1,
        )
        try:
            return (seq[_ - sorted_pos[0]].upper() for _ in position)
        except IndexError as exc:
            raise PastContigEndError(contig, max(position)) from exc
