import unittest
from reditools.region import Region
from reditools.region_collection import RegionCollection


class TestRegionCollection(unittest.TestCase):

    def setUp(self):
        self.rc = RegionCollection()
        self.rc.add_region(Region('chr1', 0, 99))
        self.rc.add_region(Region('chr1', 100, 199))
        self.rc.add_region(Region('chr2', 50, 150))

    def test_add_region_and_contains(self):
        # RegionCollection contains method requires ordered queries.
        self.assertTrue(self.rc.contains('chr1', 50))
        self.assertTrue(self.rc.contains('chr1', 150))
        self.assertFalse(self.rc.contains('chr1', 200))
        self.assertTrue(self.rc.contains('chr2', 100))
        self.assertFalse(self.rc.contains('chr2', 200))
        self.assertFalse(self.rc.contains('chrX', 1))

    def test_add_regions(self):
        regions = [
            Region('chr3', 0, 10),
            Region('chr3', 11, 20),
            Region('chr1', 200, 299),
        ]
        self.rc.add_regions(regions)
        self.assertTrue(self.rc.contains('chr3', 5))
        self.assertTrue(self.rc.contains('chr3', 15))
        self.assertFalse(self.rc.contains('chr3', 21))
        self.assertTrue(self.rc.contains('chr1', 250))
