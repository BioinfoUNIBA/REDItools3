"""Wrappers for PysamFastaFile."""

from pysam.libcfaidx import FastaFile as PysamFastaFile


class RTFastaFile(PysamFastaFile):
    """Wrapper for pysam.FastaFile to provide sequence cache."""

    def __new__(cls, *args, **kwargs):
        """
        Create a wrapper for pysam.FastaFile.

        Parameters:
            *args (list): positional arguments for PysamFastaFile constructor
            **kwargs (dict): named arguments for PysamFastaFile constructor

        Returns:
            PysamFastaFIle
        """
        return PysamFastaFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        """
        Create a wrapper for pysam.FastaFile.

        Parameters:
            *args (list): positional arguments for PysamFastaFile constructor
            **kwargs (dict): named arguments for PysamFastaFile constructor
        """
        PysamFastaFile.__init__(self)

    def get_base(self, contig, *position):
        """
        Retrieve the base at the given position.

        Parameters:
            contig (string): Chromsome name
            position (int): Zero-indexed position on reference

        Returns:
            Base the position as a string.

        Raises:
            IndexError: The position is not within the contig
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
        try:
            seq = self.fetch(contig, sorted_pos[0], sorted_pos[-1] + 1)
            return (seq[_ - sorted_pos[0]].upper() for _ in position)
        except IndexError as exc:
            raise IndexError(
                f'Base position {position} is outside the bounds of ' +
                '{contig}. Are you using the correct reference?',
            ) from exc
