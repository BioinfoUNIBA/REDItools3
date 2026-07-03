
import csv
from typing import IO, Iterator

from reditools.file_utils import open_stream
from reditools.region import Region


def _read_splice_sites(  # noqa: WPS231
        stream: IO,
) -> Iterator[tuple[str, int, str, str]]:
    reader = csv.reader(stream, delimiter=' ')
    for idx, row in enumerate(reader, start=1):
        if row[0].startswith('#'):
            continue
        try:  # noqa: WPS229
            assert len(row) == 5
            assert row[3] in ('A', 'D')
            assert row[4] in ('+', '-')
            position = int(row[1])
            yield (row[0], position, row[3], row[4])
        except (AssertionError, ValueError) as exc:
            raise ValueError(
                f'Cannot parse splice file entry ({stream.name}:{idx})'
            ) from exc

def _splice_site_to_region(
        contig: str,
        position: int,
        splice: str,
        strand: str,
        splicing_span: int,
) -> Region | None:
    strand_map = {'-': 'D', '+': 'A'}
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
    """
    Load genomic regions around splice sites from a file.

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
