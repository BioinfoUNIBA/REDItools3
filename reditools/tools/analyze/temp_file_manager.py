"""Manages temporary output files for the analyze tool."""
from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import json_args

if TYPE_CHECKING:
    from types import TracebackType
    from typing import Iterator

save_file = "region_file_list.csv"

class TempFileManager:
    """Manages the temporary output files for REDItools."""

    def __init__(self, dirpath: str, regions: list[Region] | None=None) -> None:
        """Create a new TempFileManager.

        Parameters
        ----------
        dirpath : str
            Where to save files.
        regions : list[Region]
            Analysis windows/chromosomes that need temporary files.
        """
        self.dirpath = dirpath
        if regions:
            temp_files = []
            for _ in regions:
                with tempfile.NamedTemporaryFile(
                    dir=self.dirpath,
                    delete=False,
                ) as tf:
                    temp_files.append(tf.name)
            self.region_file_list = list(zip(regions, temp_files))
            self.save_to_file()
        else:
            with Path(self.dirpath, save_file).open("r") as stream:
                reader = csv.DictReader(stream)
                self.region_file_list = [
                    (
                        Region.from_string(row["Region"]),
                        str(Path(self.dirpath, row["Filename"])),
                    )
                    for row in reader
                ]

    def __enter__(self) -> TempFileManager:
        """Open TempFileManager."""
        return self

    def __exit__(
        self,
        typ: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """If no error occurred, remove all temporary files.

        Deletes the *.done files, cli JSON, and region CSV files.
        """
        if typ is not None:
            return

        for _, filename in self.region_file_list:
            Path(f"{filename}.done").unlink()

        for temp_file in (json_args.json_args_filename, save_file):
            Path(self.dirpath, temp_file).unlink()
        try:
            Path(self.dirpath).rmdir()
        except OSError as os_exc:
            sys.stderr.write(
                "[WARNING] Could not delete temporary files directory "
                f"{self.dirpath}. {os_exc}\n",
            )

    def __iter__(self) -> Iterator[tuple[Region, str]]:
        """Iterate over region-filename pairs.

        Yields
        ------
        tuple[Region, str]
            Genomic region and path to save results to.
        """
        yield from self.region_file_list

    def __len__(self) -> int:
        """List the number of temporary output files."""
        return len(self.region_file_list)

    def concat(self, filepath: str, mode: str="w") -> None:
        """Concat all temporary files, then delete them.

        Parameters
        ----------
        filepath : str
            File to save concatenation to.
        mode : str
            Write (w) or append (a)
        """
        concat_output(
            [_[1] for _ in self.region_file_list],
            filepath,
            mode,
        )

    def save_to_file(self) -> None:
        """Save list of region files to CSV."""
        with Path(self.dirpath, save_file).open("w") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Region", "Filename"])
            for region, filename in self.region_file_list:
                writer.writerow([
                    region,
                    Path(filename).name,
                ])
