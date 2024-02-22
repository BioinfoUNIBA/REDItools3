"""Wrappers for PysamFastaFile."""

from pysam.libcfaidx import FastaFile as PysamFastaFile


class RTFastaFile(PysamFastaFile):
    """Wrapper for pysam.FastaFile to provide sequence cache."""

    def __new__(cls, *args, **kwargs):
        r"""
        Create a wrapper for pysam.FastaFile.

        Parameters:
            *args (list): positional arguments for PysamFastaFile constructor
            **kwargs (dict): named arguments for PysamFastaFile constructor

        Returns:
            PysamFastaFIle
        """
        return PysamFastaFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        r"""
        Create a wrapper for pysam.FastaFile.

        Parameters:
            *args (list): positional arguments for PysamFastaFile constructor
            **kwargs (dict): named arguments for PysamFastaFile constructor
        """
        PysamFastaFile.__init__(self)

        self._contig_name = False
        self._contig_cache = None

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
        if contig != self._contig_name:
            self._update_contig_cache(contig)
        try:
            if len(position) == 1:
                return self._contig_cache[position[0]]
            return [self._contig_cache[idx] for idx in position]
        except IndexError as exc:
            raise IndexError(
                f'Base position {position} is outside the bounds of ' +
                '{contig}. Are you using the correct reference?',
            ) from exc

    def _update_contig_cache(self, contig):
        keys = (contig, f'chr{contig}', contig.replace('chr', ''))
        for ref in keys:
            if ref in self:
                self._contig_cache = self.fetch(ref).upper()
                self._contig_name = contig
                return
        raise KeyError(f'Reference name {contig} not found in FASTA file.')
