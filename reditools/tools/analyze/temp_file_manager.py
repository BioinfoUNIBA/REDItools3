from __future__ import annotations

import csv
import os
import sys
import tempfile
from types import TracebackType
from typing import Iterator

from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import json_args

save_file = 'region_file_list.csv'

class TempFileManager:
    def __init__(self, dirpath: str, regions: list[Region] | None=None) -> None:
        self.dirpath = dirpath
        if regions:
            temp_files = [
                tempfile.NamedTemporaryFile(dir=self.dirpath, delete=False).name
                for _ in regions
            ]
            self.region_file_list = list(zip(regions, temp_files))
            with open(os.path.join(
                self.dirpath,
                save_file,
            ), 'w') as stream:
                writer = csv.writer(stream)
                writer.writerow(['Region', 'Filename'])
                for region, filename in self.region_file_list:
                    writer.writerow([region, os.path.basename(filename)])
        else:
            with open(os.path.join(self.dirpath, save_file), 'r') as stream:
                reader = csv.DictReader(stream)
                self.region_file_list = [
                    (
                        Region.from_string(row['Region']),
                        os.path.join(self.dirpath, row['Filename']),
                    )
                    for row in reader
                ]

    def __enter__(self) -> TempFileManager:
        return self

    def __iter__(self) -> Iterator:
        yield from self.region_file_list

    def __len__(self) -> int:
        return len(self.region_file_list)

    def concat(self, filepath: str, mode: str='w') -> None:
        concat_output(
            [_[1] for _ in self.region_file_list],
            filepath,
            mode,
        )

    def __exit__(
        self,
        exc_type: type,
        exc_value: Exception,
        traceback: TracebackType,
    ) -> None:
        if exc_type is not None:
            return

        for _, filename in self.region_file_list:
            os.remove(f'{filename}.done')

        for temp_file in (json_args.json_args_filename, save_file):
            os.remove(os.path.join(self.dirpath, temp_file))
        try:
            os.rmdir(self.dirpath)
        except OSError as exc:
            sys.stderr.write(
                '[WARNING] Could not delete temporary files directory '
                f'{self.dirpath}. {exc}\n'
            )
