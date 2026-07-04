import os
import unittest
from tempfile import NamedTemporaryFile

from reditools.region import Region
from reditools.splicing_file import load_splicing_file


class TestSplicingFile(unittest.TestCase):
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
        splice_sites = list(load_splicing_file(fname, 5))
        self.check_test_data(test_data, splice_sites)
        os.remove(fname)

    def test_splicing_edge(self):
        test_data = [
            ('chr1', '1', '25', 'A', '+'),
            ('chr1', '1', '25', 'D', '-'),
            ('chr1', '3', '25', 'D', '-'),
        ]
        fname = self.write_file(['#Header'] + test_data)
        splice_sites = list(load_splicing_file(fname, 5))
        self.assertEqual(
            splice_sites,
            [Region(contig='chr1', start=0, stop=2)],
        )
        os.remove(fname)
