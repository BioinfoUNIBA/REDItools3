
import csv
import os
from gzip import open as gzip_open
from typing import IO, Iterator

from reditools.region import Region


def open_stream(path: str, mode: str='rt', encoding: str='utf-8'):
    """
    Open a file stream, handling both plain and gzipped files.

    Parameters
    ----------
    path : str
        The path to the file.
    mode : str, optional
        The mode in which the file is opened (default is 'rt').
    encoding : str, optional
        The encoding to use (default is 'utf-8').

    Returns
    -------
    file-like object
        The opened file stream.
    """
    if path.endswith('gz'):
        return gzip_open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)  # noqa: WPS515


def read_bed_file(*path: str) -> Iterator[Region]:
    """
    Read genomic regions from one or more BED files.

    Parameters
    ----------
    *path : str
        Paths to the BED files.

    Yields
    ------
    Region
        The regions defined in the BED files.
    """
    if len(path) > 1:
        yield from read_bed_file(*path[1:])
    with open_stream(path[0]) as stream:
        reader = csv.reader(
            filter(lambda row: row[0] != '#', stream),
            delimiter='\t',
        )
        for row in reader:
            yield Region(
                contig=row[0],
                start=int(row[1]),
                stop=int(row[2]),
            )


def concat(
        output: IO,
        *fnames: str,
        clean_up: bool=True,
        encoding: str='utf-8',
) -> None:
    """
    Concatenate multiple files into a single output stream.

    Parameters
    ----------
    output : IO
        The output stream to write to.
    *fnames : str
        The names of the files to concatenate.
    clean_up : bool, optional
        Whether to remove the source files after concatenation
        (default is True).
    encoding : str, optional
        The encoding to use when reading files (default is 'utf-8').
    """
    for fname in fnames:
        with open(fname, 'r', encoding=encoding) as stream:
            for line in stream:
                output.write(line)
        if clean_up:
            os.remove(fname)


def load_text_file(file_name: str) -> list[str]:
    """
    Load lines from a text file into a list, stripping whitespace.

    Parameters
    ----------
    file_name : str
        The name of the file to load.

    Returns
    -------
    list[str]
        A list of stripped lines from the file.
    """
    with open_stream(file_name, 'r') as stream:
        return [line.strip() for line in stream]


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
