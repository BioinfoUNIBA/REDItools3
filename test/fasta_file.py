import unittest
from tempfile import NamedTemporaryFile
import os
from reditools.fasta_file import RTFastaFile
import random
from itertools import chain


class TestRTFastaFile(unittest.TestCase):
    def setUp(self):
        self.contig1 = 'test1'
        self.seq1 = ''.join([random.choice('ACTG') for _ in range(80)])
        self.contig2 = 'chrtest2'
        self.seq2 = ''.join([random.choice('ACTG') for _ in range(80)])
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8') as f:
            self.fasta_fname = f.name
            f.write(f'>{self.contig1}\n{self.seq1}\n')
            f.write(f'>{self.contig2}\n{self.seq2}\n')

    def tearDown(self):
        os.remove(self.fasta_fname)

    def test_get_base(self):
        with RTFastaFile(self.fasta_fname) as rff:
            positions = list(range(len(self.seq1)))
            fasta_seq = rff.get_base(self.contig1, *positions)
            self.assertEqual(self.seq1, ''.join(fasta_seq))

    def test_get_base_splice(self):
        with RTFastaFile(self.fasta_fname) as rff:
            positions = list(chain(
                range(20),
                range(len(self.seq1) - 20, len(self.seq1)),
            ))
            fasta_seq = rff.get_base(self.contig1, *positions)
            self.assertEqual(
                self.seq1[:20] + self.seq1[-20:],
                ''.join(fasta_seq),
            )

    def test_get_base_prefix(self):
        with RTFastaFile(self.fasta_fname) as rff:
            positions = range(len(self.seq2))
            fasta_seq = rff.get_base('test2', *positions)
            self.assertEqual(self.seq2, ''.join(fasta_seq))

            fasta_seq = rff.get_base('chrtest2', *positions)
            self.assertEqual(self.seq2, ''.join(fasta_seq))

    def test_get_base_missing_contig(self):
        with RTFastaFile(self.fasta_fname) as rff:
            with self.assertRaises(KeyError):
                rff.get_base('test3', 0)
            with self.assertRaises(KeyError):
                rff.get_base('chrtest3', 0)

    def test_get_base_out_of_bounds(self):
        with RTFastaFile(self.fasta_fname) as rff:
            with self.assertRaises(IndexError):
                positions = range(len(self.seq1) - 20, len(self.seq1) + 20)
                seq_iter = rff.get_base(
                    self.contig1,
                    *positions,
                )
                list(seq_iter)
