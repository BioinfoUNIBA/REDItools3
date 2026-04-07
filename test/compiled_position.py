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

    def test_alts_and_variants(self):
        self.cp.add_base(10, '+', 'C')
        self.cp.add_base(10, '+', 'A')
        self.assertEqual(self.cp.alts, ['C'])
        self.assertEqual(self.cp.variants, ['AC'])

    def test_calculate_strand(self):
        self.cp.add_base(20, '+', 'A')
        self.cp.add_base(21, '-', 'A')
        self.cp.add_base(22, '+', 'C')
        with self.assertRaises(ValueError):
            self.cp.strand
        self.assertEqual(self.cp.calculate_strand(), '+')
        self.assertEqual(self.cp.strand, '+')
        self.assertEqual(self.cp.calculate_strand(0.7), '*')
        self.assertEqual(self.cp.strand, '*')

    def test_filter_by_strand(self):
        self.cp.add_base(5, '+', 'A')
        self.cp.add_base(5, '+', 'A')
        self.cp.add_base(6, '-', 'C')
        self.cp.calculate_strand()
        self.cp.filter_by_strand()
        self.assertEqual(len(self.cp), 2)
        self.assertEqual(self.cp['A'], 2)
        self.assertEqual(self.cp['C'], 0)
        self.assertEqual(self.cp.strand, '+')

    def test_filter_by_strand_star(self):
        self.cp.add_base(5, '*', 'A')
        self.cp.add_base(5, '*', 'A')
        self.cp.add_base(6, '+', 'C')
        self.cp.add_base(6, '-', 'T')
        self.assertEqual(
            self.cp.calculate_strand(threshold=1),
            '*',
        )
        self.cp.filter_by_strand()
        self.assertEqual(len(self.cp), 4)
        self.assertEqual(self.cp['A'], 2)
        self.assertEqual(self.cp['C'], 1)

    def test_reference(self):
        self.assertEqual(self.cp.reference, 'A')

    def test_mean_quality(self):
        self.assertEqual(self.cp.mean_quality, 0)
        self.cp.add_base(5, '*', 'A')
        self.cp.add_base(6, '+', 'C')
        self.assertEqual(self.cp.mean_quality, 5.5)

    def test_edit_ratio(self):
        self.assertEqual(self.cp.edit_ratio, 0)
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'C')
        self.cp.add_base(30, '+', 'C')
        self.assertEqual(self.cp.edit_ratio, 0.75)

    def test_depth(self):
        self.assertEqual(self.cp.depth, 0)
        self.assertEqual(self.cp.depth, len(self.cp))
        self.cp.add_base(40, '+', 'A')
        self.assertEqual(self.cp.depth, 1)
        self.assertEqual(self.cp.depth, len(self.cp))

    def test_per_base_depth(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'C')
        self.assertEqual(self.cp.per_base_depth, [1, 2, 0, 0])
