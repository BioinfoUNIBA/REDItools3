"""Class to annotate RNA editing sites with DNA data."""
from __future__ import annotations

import csv
from typing import IO, Iterator

from reditools import file_utils
from reditools.compiled_position import RTResult
from reditools.constants import comp_map


class AnalyzeMismatchError(ValueError):
    """Reference bases from two REDItools output files do not match."""

    def __init__(self) -> None:
        """Initialize self."""
        self.message = "Files do not appear to use the same reference."
        super().__init__(self.message)

def read_rta_file(file_stream: IO) -> Iterator[RTResult]:
    """Read an Analayze tool output.

    Parameters
    ----------
    file_stream : IO
        Stream to read from.

    Yields
    ------
    RTResult
    """
    reader = csv.DictReader(file_stream, delimiter="\t")
    for row in reader:
        yield RTResult.from_dict(row)

class RTAnnotater:
    """Class to annotate RNA editing sites with DNA data.

    It merges RNA and DNA editing tables by matching positions.
    """

    def __init__(
        self,
        contig_order: dict[str, int],
        do_complement: bool=False,
    ) -> None:
        """Initialize RTAnnotater.

        Parameters
        ----------
        contig_order : dict[str, int]
            A dictionary mapping contig names to their sort order.
        do_complement : bool
            If True, annotate() will report the DNA complement if the RNA data
            comes from the minus strand.
        """
        self.contig_order = contig_order
        self.do_complement = do_complement

    def annotate(
        self,
        rna_iter: Iterator[RTResult],
        dna_iter: Iterator[RTResult],
    ) -> Iterator[RTResult]:
        """Read input files and yields annotated results.

        Parameters
        ----------
        rna_iter : Iterator[RTResult]
            Iterator of RTResults from RNA analysis.
        dna_iter : Iterator[RTResult]
            Iterator of RTResults from analysis to use as annotation.

        Yields
        ------
        RTResult
            rna_iter results annotated with dna_iter results.
        """
        dna_entry = next(dna_iter, None)

        for rna_entry in rna_iter:
            while dna_entry is not None and \
                    self.cmp_position(rna_entry, dna_entry) > 0:
                dna_entry = next(dna_iter, None)
            if dna_entry is not None and \
                    self.cmp_position(rna_entry, dna_entry) == 0:
                self.annotate_row(rna_entry, dna_entry)
            yield rna_entry

    def filtered_annotate(
        self,
        rna_iter: Iterator[RTResult],
        dna_iter: Iterator[RTResult],
    ) -> Iterator[RTResult]:
        """Read input files and yields annotated results.

        Works the same as annotate(), except only entries both DNA and RNA
        reads are yielded..

        Parameters
        ----------
        rna_iter : Iterator[RTResult]
            Iterator of RTResults from RNA analysis.
        dna_iter : Iterator[RTResult]
            Iterator of RTResults from analysis to use as annotation.

        Yields
        ------
        RTResult
            rna_iter results annotated with dna_iter results.
        """
        dna_entry = next(dna_iter, None)

        for rna_entry in rna_iter:
            while dna_entry is not None and \
                    self.cmp_position(rna_entry, dna_entry) > 0:
                dna_entry = next(dna_iter, None)
            if dna_entry is not None and \
                    self.cmp_position(rna_entry, dna_entry) == 0:
                self.annotate_row(rna_entry, dna_entry)
                yield rna_entry

    def cmp_position(
            self,
            rna_entry: RTResult,
            dna_entry: RTResult,
    ) -> int:
        """Compare the positions of RNA and DNA entries.

        Parameters
        ----------
        rna_entry : RTResult
            A row from the RNA editing file.
        dna_entry : RTResult
            A row from the DNA editing file, or None if the end is reached.

        Returns
        -------
        int
            Negative if RNA < DNA, positive if RNA > DNA, zero if equal.
        """
        if rna_entry.contig == dna_entry.contig:
            return rna_entry.position - dna_entry.position
        rna_contig_idx = self.contig_order[rna_entry.contig]
        dna_contig_idx = self.contig_order[dna_entry.contig]
        return rna_contig_idx - dna_contig_idx

    def annotate_row(
            self,
            rna_row: RTResult,
            dna_row: RTResult,
    ) -> None:
        """Add DNA information to an RNA row.

        Parameters
        ----------
        rna_row : RTResult
            A row from the RNA editing file.
        dna_row : RTResult
            A matching row from the DNA editing file.
        """
        if rna_row.reference == comp_map[dna_row.reference]:
            if self.do_complement:
                dna_row.complement()
        elif rna_row.reference != dna_row.reference:
            raise AnalyzeMismatchError
        rna_row.dna = dna_row

    def merge_files(
            self,
            rna_file: str,
            dna_file: str,
    ) -> Iterator[RTResult]:
        """Merge RNA and DNA files and yield annotated rows.

        Parameters
        ----------
        rna_file : str
            Path to the RNA editing file.
        dna_file : str
            Path to the DNA editing file.

        Yields
        ------
        RTResult
            Annotated (or original if no match) RNA row.
        """
        with file_utils.open_stream(rna_file, "r") as rna_stream, \
                file_utils.open_stream(dna_file, "r") as dna_stream:
            rna_reader = read_rta_file(rna_stream)
            dna_reader = read_rta_file(dna_stream)
            yield from self.annotate(rna_reader, dna_reader)
