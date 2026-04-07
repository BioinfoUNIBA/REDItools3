import unittest
from reditools.compiled_position import CompiledPosition
from reditools.rtresult import RTResult


class TestRTResult(unittest.TestCase):
    def setUp(self):
        self.cp = CompiledPosition('A', 'chr1', 100)

    def test_base_properties(self):
        result = RTResult(self.cp, '+')
        self.assertEqual(result.contig, 'chr1')
        self.assertEqual(result.position, 101)
        self.assertEqual(result.reference, 'A')

    def test_no_variants(self):
        result = RTResult(self.cp, '+')
        self.assertEqual(result.variants, [])

    def test_variants(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'G')
        result = RTResult(self.cp, '+')
        self.assertEqual(sorted(result.variants), ['AC', 'AG'])

    def test_mean_quality(self):
        result = RTResult(self.cp, '+')
        self.assertEqual(result.mean_quality, 0)

        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'G')
        result = RTResult(self.cp, '+')
        self.assertEqual(result.mean_quality, 35)

    def test_edit_ratio(self):
        result = RTResult(self.cp, '+')
        self.assertEqual(result.edit_ratio, 0)

        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'C')
        self.cp.add_base(30, '+', 'C')
        result = RTResult(self.cp, '+')
        self.assertEqual(result.edit_ratio, 0.75)

    def test_depth(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        result = RTResult(self.cp, '+')
        self.assertEqual(result.depth, 2)

    def test_per_base_depth(self):
        self.cp.add_base(40, '+', 'A')
        self.cp.add_base(35, '-', 'C')
        self.cp.add_base(30, '+', 'C')
        result = RTResult(self.cp, '+')
        self.assertEqual(result.per_base_depth, [1, 2, 0, 0])
