import unittest
from reditools.compiled_position import CompiledPosition


class TestCompiledPosition(unittest.TestCase):
    def setUp(self):
        self.cp = CompiledPosition('A', 'chr1', 100)

    def test_add_base_and_len(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'G')
        self.assertEqual(len(self.cp), 3)

    def test_get_base_counts(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'A')
        self.cp.add_base(30, '+', 'C')
        self.assertEqual(self.cp['A'], 2)
        self.assertEqual(self.cp['C'], 1)
        self.assertEqual(self.cp['REF'], 2)

    def test_iter(self):
        self.cp.add_base(41, '+', 'A')
        self.cp.add_base(42, '+', 'C')
        self.cp.add_base(43, '+', 'G')
        counts = list(self.cp)
        self.assertEqual(counts, [1, 1, 1, 0])

    def test_complement(self):
        self.cp.add_base(11, '+', 'A')
        self.cp.add_base(12, '-', 'C')
        self.cp.complement()
        self.assertEqual(self.cp.bases, ['T', 'G'])
        self.assertEqual(self.cp.ref, 'T')

    def test_get_variants(self):
        self.cp.add_base(10, '+', 'C')
        self.cp.add_base(10, '+', 'A')
        variants = self.cp.get_variants()
        self.assertEqual(variants, ['C'])

    def test_get_strand(self):
        self.cp.add_base(20, '+', 'A')
        self.cp.add_base(21, '-', 'A')
        self.cp.add_base(22, '+', 'C')
        self.assertEqual(self.cp.get_strand(), '+')
        self.assertEqual(self.cp.get_strand(0.7), '*')

    def test_filter_by_strand(self):
        self.cp.add_base(5, '+', 'A')
        self.cp.add_base(6, '-', 'C')
        self.cp.filter_by_strand('+')
        self.assertEqual(len(self.cp), 1)
        self.assertEqual(self.cp['A'], 1)
        self.assertEqual(self.cp['C'], 0)
        self.assertEqual(self.cp.get_strand(1), '+')

    def test_filter_by_strand_star(self):
        self.cp.add_base(5, '*', 'A')
        self.cp.add_base(6, '+', 'C')
        self.cp.filter_by_strand('*')
        self.assertEqual(len(self.cp), 2)
        self.assertEqual(self.cp['A'], 1)
        self.assertEqual(self.cp['C'], 1)
