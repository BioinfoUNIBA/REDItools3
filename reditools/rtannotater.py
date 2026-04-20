import csv
from typing import IO, Any, Iterator

from reditools import file_utils


class RTAnnotater:
    """Class to annotate RNA editing sites with DNA data.

    It merges RNA and DNA editing tables by matching positions.

    Attributes
    ----------
    legacy_map : tuple
        Mapping of old field names to new field names for backward
        compatibility.
    """

    legacy_map = (
        ('Coverage-q30', 'Coverage'),
        ('gCoverage-q30', 'gCoverage'),
    )

    def __init__(
            self,
            rna_file: str,
            dna_file: str,
            contig_order: dict[str, int],
    ):
        """Initialize RTAnnotater.

        Parameters
        ----------
        rna_file : str
            Path to the RNA editing file.
        dna_file : str
            Path to the DNA editing file.
        contig_order : dict[str, int]
            A dictionary mapping contig names to their sort order.
        """
        self.rna_file = rna_file
        self.dna_file = dna_file
        self.contig_order = contig_order

    def annotate(self, stream: IO):
        """Read input files and write annotated results to a stream.

        Parameters
        ----------
        stream : IO
            The output stream to write the annotated table.
        """
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=[
            'Region',
            'Position',
            'Reference',
            'Strand',
            'Coverage',
            'MeanQ',
            'BaseCount[A,C,G,T]',
            'AllSubs',
            'Frequency',
            'gCoverage',
            'gMeanQ',
            'gBaseCount[A,C,G,T]',
            'gAllSubs',
            'gFrequency'])
        writer.writeheader()
        writer.writerows(self.merge_files())

    def cmp_position(
            self,
            rna_entry: dict[Any, Any],
            dna_entry: dict[Any, Any] | None,
    ) -> int:
        """Compare the positions of RNA and DNA entries.

        Parameters
        ----------
        rna_entry : dict[Any, Any]
            A row from the RNA editing file.
        dna_entry : dict[Any, Any] | None
            A row from the DNA editing file, or None if the end is reached.

        Returns
        -------
        int
            Negative if RNA < DNA, positive if RNA > DNA, zero if equal.
        """
        if dna_entry is None:
            return -1
        rna_contig_idx = self.contig_order[rna_entry['Region']]
        # If the DNA contig is not in the RNA file, assume its position is
        # earlier than the current RNA contig to induce fast-forwarding.
        dna_contig_idx = self.contig_order.get(
            dna_entry['Region'],
            0
        )
        if rna_contig_idx == dna_contig_idx:
            return int(rna_entry['Position']) - int(dna_entry['Position'])
        return rna_contig_idx - dna_contig_idx

    def annotate_row(
            self,
            rna_row: dict[Any, Any],
            dna_row: dict[Any, Any],
    ) -> dict[Any, Any]:
        """Add DNA information to an RNA row.

        Parameters
        ----------
        rna_row : dict[Any, Any]
            A row from the RNA editing file.
        dna_row : dict[Any, Any]
            A matching row from the DNA editing file.

        Returns
        -------
        dict[Any, Any]
            The annotated RNA row.
        """
        rna_row['gCoverage'] = dna_row['Coverage']
        rna_row['gMeanQ'] = dna_row['MeanQ']
        rna_row['gBaseCount[A,C,G,T]'] = dna_row['BaseCount[A,C,G,T]']
        rna_row['gAllSubs'] = dna_row['AllSubs']
        rna_row['gFrequency'] = dna_row['Frequency']
        return rna_row

    def legacy_translate(self, row: dict[Any, Any]) -> dict[Any, Any]:
        """Translate legacy field names to current ones.

        Parameters
        ----------
        row : dict[Any, Any]
            A row from an editing file.

        Returns
        -------
        dict[Any, Any]
            The translated row.
        """
        for old_key, new_key in self.legacy_map:
            if old_key in row:
                row[new_key] = row.pop(old_key)
        return row

    def merge_files(self) -> Iterator[dict[Any, Any]]:
        """Merge RNA and DNA files and yield annotated rows.

        Yields
        ------
        dict[Any, Any]
            Annotated (or original if no match) RNA row.
        """
        with file_utils.open_stream(self.rna_file, 'r') as rna_stream, \
                file_utils.open_stream(self.dna_file, 'r') as dna_stream:
            rna_reader = csv.DictReader(rna_stream, delimiter='\t')
            dna_reader = csv.DictReader(dna_stream, delimiter='\t')

            dna_entry = next(dna_reader, None)

            for rna_entry in rna_reader:
                self.legacy_translate(rna_entry)

                while self.cmp_position(rna_entry, dna_entry) > 0:
                    dna_entry = next(dna_reader, None)
                if dna_entry is not None and \
                        self.cmp_position(rna_entry, dna_entry) == 0:
                    self.legacy_translate(dna_entry)
                    yield self.annotate_row(rna_entry, dna_entry)
                else:
                    yield rna_entry
