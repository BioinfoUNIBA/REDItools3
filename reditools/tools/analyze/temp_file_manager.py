from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Iterator

from reditools.region import Region
from reditools.tools.analyze.concat_output import concat_output
from reditools.tools.analyze.parse_args import json_args

json_windows_file = 'tempfile_map.json'

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
                json_windows_file,
            ), 'w') as stream:
                json.dump(self.region_file_list, stream, default=str)
        else:
            self._load_from_file()

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

    def cleanup(self) -> None:
        for _, filename in self.region_file_list:
            os.remove(f'{filename}.done')

        for temp_file in (json_args.json_args_filename, json_windows_file):
            os.remove(os.path.join(self.dirpath, temp_file))
        try:
            os.rmdir(self.dirpath)
        except OSError as exc:
            sys.stderr.write(
                '[WARNING] Could not delete temporary files directory '
                f'{self.dirpath}. {exc}\n'
            )

    def _load_from_file(self) -> None:
        with open(os.path.join(self.dirpath, json_windows_file), 'r') as stream:
            self.region_file_list = json.load(stream)
