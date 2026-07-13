"""A wrapper around pysam.FastaFile for genomic sequence access."""
from __future__ import annotations

from typing import TYPE_CHECKING

from pysam.libcfaidx import FastaFile as PysamFastaFile

if TYPE_CHECKING:
    from types import TracebackType
    from typing import Iterator


class MissingContigError(LookupError):
    """Contig name is missing from the FASTA file."""

    def __init__(self, contig_name: str) -> None:
        """Initialize self.

        Parameters
        ----------
        contig_name : str
            Missing contig name.
        """
        self.message = f"Reference name {contig_name} not found in FASTA file."
        super().__init__(self.message)

class PastContigEndError(LookupError):
    """Genomic position is outside contig bounds."""

    def __init__(self, contig_name: str, position: int) -> None:
        """Initialize self.

        Parameters
        ----------
        contig_name : str
            Name of the contig.
        position : int
            Offending genomic position.
        """
        self.message = (
            f"Base position {position} is outside the bounds of "
            f"{contig_name}. Are you using the correct reference?"
        )
        super().__init__(self.message)

class RTFastaFile:
    """A wrapper around pysam.FastaFile for genomic sequence access."""

    def __init__(self, filename: str) -> None:
        """Initialize the RTFastaFile.

        Parameters
        ----------
        filename : str
            FASTA file path.
        """
        self.pysam_fasta_file = PysamFastaFile(filename)

    def __enter__(self) -> RTFastaFile:
        """Open RTFastaFile."""
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
        if contig not in self.pysam_fasta_file.references:
            if contig.startswith("chr"):
                new_contig = contig[3:]
            else:
                new_contig = f"chr{contig}"
            if new_contig not in self.pysam_fasta_file.references:
                raise MissingContigError(contig)
            contig = new_contig
        sorted_pos = sorted(position)
        seq = self.pysam_fasta_file.fetch(
            contig,
            sorted_pos[0],
            sorted_pos[-1] + 1,
        )
        try:
            for pos in position:
                yield seq[pos - sorted_pos[0]].upper()
        except IndexError as exc:
            raise PastContigEndError(contig, max(position)) from exc
