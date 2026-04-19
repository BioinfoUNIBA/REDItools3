from typing import Iterator

from pysam.libcfaidx import FastaFile as PysamFastaFile


class RTFastaFile(PysamFastaFile):
    """
    A wrapper around pysam.FastaFile for genomic sequence access.
    """

    def __new__(cls, *args, **kwargs):
        """
        Create a new instance of RTFastaFile.

        Parameters
        ----------
        *args
            Arguments passed to pysam.FastaFile.
        **kwargs
            Keyword arguments passed to pysam.FastaFile.

        Returns
        -------
        RTFastaFile
            A new instance of RTFastaFile.
        """
        return PysamFastaFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        """
        Initialize the RTFastaFile.

        Parameters
        ----------
        *args
            Arguments passed to pysam.FastaFile.
        **kwargs
            Keyword arguments passed to pysam.FastaFile.
        """
        PysamFastaFile.__init__(self)

    def get_base(self, contig: str, *position: int) -> Iterator[str]:
        """
        Retrieve bases at specified positions from a contig.

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
        KeyError
            If the contig is not found in the FASTA file.
        IndexError
            If a position is outside the bounds of the contig.
        """

        if contig not in self:
            if contig.startswith('chr'):
                new_contig = contig.replace('chr', '')
            else:
                new_contig = f'chr{contig}'
            if new_contig not in self:
                raise KeyError(
                    f'Reference name {contig} not found in FASTA file.',
                )
            contig = new_contig
        sorted_pos = sorted(position)
        seq = self.fetch(
            contig,
            sorted_pos[0],
            sorted_pos[-1] + 1,
        )
        try:
            return (seq[_ - sorted_pos[0]].upper() for _ in position)
        except IndexError as exc:
            raise IndexError(
                f'Base position {position} is outside the bounds of ' +
                '{contig}. Are you using the correct reference?',
            ) from exc
