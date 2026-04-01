import unittest
from tempfile import NamedTemporaryFile
import os
import reditools.file_utils as file_utils
from reditools.region import Region
import gzip


class TestFileUtils(unittest.TestCase):
    def test_open_stream_plain(self):
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            f.write('test123')
            fname = f.name
        with file_utils.open_stream(fname, 'rt') as stream:
            content = stream.read()
        self.assertEqual(content, 'test123')
        os.remove(fname)

    def test_open_stream_gzip(self):
        with NamedTemporaryFile(delete=False,
                                suffix='.gz',
                                mode='wb') as f:
            f.write(gzip.compress(b'test_gzip'))
            fname = f.name
        with file_utils.open_stream(fname, 'rt') as stream:
            content = stream.read()
        self.assertEqual('test_gzip', content)
        os.remove(fname)

    def test_read_bed_file(self):
        bed_lines = ["chr1\t10\t20", "chr2\t30\t40"]
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            fname = f.name
            for line in bed_lines:
                f.write(line + "\n")
        f.close()

        result = list(file_utils.read_bed_file(fname))
        self.assertEqual(result[0], Region('chr1', 10, 20))
        self.assertEqual(result[1], Region('chr2', 30, 40))
        os.remove(fname)

    def test_concat(self):
        file_contents = ['file1', 'file2', 'file3']
        file_names = []
        for content in file_contents:
            with NamedTemporaryFile(delete=False,
                                    mode='w',
                                    encoding='utf-8') as stream:
                file_names.append(stream.name)
                stream.write(f'{content}\n')

        output = NamedTemporaryFile(delete=False, mode='w+', encoding='utf-8')
        file_utils.concat(output, *file_names, encoding='utf-8')
        output.seek(0)
        result = output.read()
        self.assertEqual(result, ''.join([f'{_}\n' for _ in file_contents]))

        for fname in file_names:
            self.assertFalse(os.path.exists(fname))
        output.close()
        os.remove(output.name)

    def test_load_text_file(self):
        text_lines = ["rowA", "rowB", "rowC"]
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            fname = f.name
            f.write('\n'.join(text_lines))
        result = file_utils.load_text_file(fname)
        self.assertEqual(result, text_lines)
        os.remove(fname)

    def test_splicing_basic(self):
        data = [
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
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            fname = f.name
            f.write('#Header\n')
            for entry, _ in data:
                f.write(' '.join(entry))
                f.write('\n')
        test_iter = zip(
            file_utils.load_splicing_file(fname, 5),
            data,
        )
        for splice_region, (_, test_region) in test_iter:
            self.assertEqual(splice_region, test_region)
        os.remove(fname)

    def test_splicing_edge(self):
        data = [
            ('chr1', '1', '25', 'A', '+'),
            ('chr1', '1', '25', 'D', '-'),
            ('chr1', '3', '25', 'D', '-'),
        ]
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            fname = f.name
            f.write('#Header\n')
            for entry in data:
                f.write(' '.join(entry))
                f.write('\n')
        splice_sites = list(file_utils.load_splicing_file(fname, 5))
        self.assertEqual(
            splice_sites,
            [Region(contig='chr1', start=0, stop=2)],
        )
        os.remove(fname)
