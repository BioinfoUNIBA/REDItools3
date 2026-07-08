"""Load genomic regions around splice sites from a file."""
from __future__ import annotations

import csv
from typing import IO, Iterator

from reditools.file_utils import open_stream
from reditools.region import Region


class SpliceFileFormatError(ValueError):
    """Splice file is not in expected format."""

    def __init__(self, file_name: str, line_number: int) -> None:
        """Initialize self.

        Parameters
        ----------
        file_name : str
            The offending file.
        line_number : int
            The line number that does not match expected format.
        """
        self.message = (
            f"Cannot parse splice file entry ({file_name}:{line_number})"
        )
        super().__init__(self.message)

def _read_splice_sites(  # noqa: WPS231
        stream: IO,
) -> Iterator[tuple[str, int, str, str]]:
    reader = csv.reader(stream, delimiter=" ")
    for idx, row in enumerate(reader, start=1):
        if row[0].startswith("#"):
            continue
        if len(row) != 5 or \
                row[3] not in ("A", "D") or \
                row[4] not in ("+", "-"):  # noqa: PLR2004
            raise SpliceFileFormatError(stream.name, idx)
        try:
            position = int(row[1])
        except ValueError as exc:
            raise SpliceFileFormatError(stream.name, idx) from exc
        yield (row[0], position, row[3], row[4])

def _splice_site_to_region(
        contig: str,
        position: int,
        splice: str,
        strand: str,
        splicing_span: int,
) -> Region | None:
    strand_map = {"-": "D", "+": "A"}
    position = position - 1
    if strand_map[strand] == splice:
        start = max(position - splicing_span, 0)
        stop = position
    else:
        start = position
        stop = position + splicing_span
    if start != stop:
        return Region(contig=contig, start=start, stop=stop)
    return None

def load_splicing_file(
        splicing_file: str,
        splicing_span: int,
) -> Iterator[Region]:
    """Load genomic regions around splice sites from a file.

    Splice site files are space delimited and have five columns:
    1. Chromosome/contig
    2. Genomic start
    3. Genomic stop (ignored)
    4. Splice type (A or D)
    5. Strand (+ or -)

    Parameters
    ----------
    splicing_file : str
        The path to the splice sites file.
    splicing_span : int
        The number of bases around each splice site to include in the region.

    Yields
    ------
    Region
        The genomic regions around the splice sites.
    """
    with open_stream(splicing_file) as stream:
        for splice_data in _read_splice_sites(stream):
            region = _splice_site_to_region(*splice_data, splicing_span)
            if region is not None:
                yield region
