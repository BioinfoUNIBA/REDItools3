import gzip
import os
import unittest
from tempfile import NamedTemporaryFile

from reditools import file_utils
from reditools.region import Region


class TestFileUtils(unittest.TestCase):
    def write_file(self, data_list, sep=' '):
        with NamedTemporaryFile(
                delete=False,
                mode='w',
                encoding='utf-8',
        ) as stream:
            for row in data_list:
                if isinstance(row, str):
                    stream.write(row)
                else:
                    stream.write(sep.join([str(_) for _ in row]))
                stream.write('\n')
            return stream.name

    def check_test_data(self, test_data, real_data):
        self.assertEqual([_[1] for _ in test_data], real_data)

    def test_open_stream_plain(self):
        test_str = 'test123'
        with NamedTemporaryFile(
                delete=False,
                mode='w',
                encoding='utf-8',
        ) as stream:
            stream.write(test_str)
            fname = stream.name
        with file_utils.open_stream(fname, 'rt') as stream:
            file_content = stream.read()
        self.assertEqual(file_content, test_str)
        os.remove(fname)

    def test_open_stream_gzip(self):
        test_str = 'test_gzip'
        with NamedTemporaryFile(
                delete=False,
                suffix='.gz',
                mode='wb',
        ) as stream:
            stream.write(gzip.compress(bytes(test_str, 'utf-8')))
            fname = stream.name
        with file_utils.open_stream(fname, 'rt') as stream:
            file_content = stream.read()
        self.assertEqual(file_content, test_str)
        os.remove(fname)

    def test_read_bed_file(self):
        bed_data = (
            (
                ('chr1', 10, 20),
                Region('chr1', 10, 20),
            ),
            (
                ('chr2', 30, 40),
                Region('chr2', 30, 40),
            ),
        )
        fname = self.write_file((_[0] for _ in bed_data), sep='\t')
        region_list = list(file_utils.read_bed_file(fname))
        self.check_test_data(bed_data, region_list)
        os.remove(fname)

    def test_concat(self):
        file_contents = ('file1', 'file2', 'file3')
        file_names = [self.write_file([_]) for _ in file_contents]

        with NamedTemporaryFile(
                delete=False,
                mode='w',
                encoding='utf-8') as stream:
            file_utils.concat(stream, *file_names, encoding='utf-8')
            concat_filename = stream.name
        for fname in file_names:
            self.assertFalse(os.path.exists(fname))

        with open(concat_filename, 'r') as stream:
            self.assertEqual(
                stream.read(),
                ''.join([f'{_}\n' for _ in file_contents]),
            )

        os.remove(concat_filename)

    def test_load_text_file(self):
        text_lines = ["rowA", "rowB", "rowC"]
        with NamedTemporaryFile(
                delete=False,
                mode='w',
                encoding='utf-8',
        ) as stream:
            fname = stream.name
            stream.write('\n'.join(text_lines))
        loaded_text = file_utils.load_text_file(fname)
        self.assertEqual(loaded_text, text_lines)
        os.remove(fname)

    def test_splicing_basic(self):
        test_data = [
            (
                ('chr1', '10', '25', 'A', '+'),
                Region(contig='chr1', start=4, stop=9),
            ),
            (
                ('chr2', '20', '25', 'D', '-'),
                Region(contig='chr2', start=14, stop=19),
            ),
            (
                ('chr3', '5', '15', 'A', '-'),
                Region(contig='chr3', start=4, stop=9),
            ),
            (
                ('chr3', '5', '10', 'D', '+'),
                Region(contig='chr3', start=4, stop=9),
            ),
        ]
        fname = self.write_file(
            ['#Header'] + [_[0] for _ in test_data],
        )
        splice_sites = list(file_utils.load_splicing_file(fname, 5))
        self.check_test_data(test_data, splice_sites)
        os.remove(fname)

    def test_splicing_edge(self):
        test_data = [
            ('chr1', '1', '25', 'A', '+'),
            ('chr1', '1', '25', 'D', '-'),
            ('chr1', '3', '25', 'D', '-'),
        ]
        fname = self.write_file(['#Header'] + test_data)
        splice_sites = list(file_utils.load_splicing_file(fname, 5))
        self.assertEqual(
            splice_sites,
            [Region(contig='chr1', start=0, stop=2)],
        )
        os.remove(fname)
