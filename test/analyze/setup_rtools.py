import unittest

from reditools.reditools import REDItools
from reditools.tools.analyze.parse_args import parse_args
from reditools.tools.analyze.setup_rtools import setup_rtools


class TestSetupRTools(unittest.TestCase):
    def test_options(self):
        options = parse_args([
            'example.bam',
            '-r', 'genome.fa',
            '--min-base-position', '10',
            '--max-base-position', '30',
            '--min-base-quality', '23',
            '--strand', '1',
            '--strand-confidence-threshold', '0.567',
            '--strand-correction',
        ])
        rtools = setup_rtools(options)
        self.assertIsInstance(rtools, REDItools)
        self.assertEqual(rtools.reference, 'genome.fa')
        self.assertEqual(rtools.min_base_position, 10)
        self.assertEqual(rtools.max_base_position, 30)
        self.assertEqual(rtools.min_base_quality, 23)
        self.assertEqual(rtools.strand, 1)
        self.assertEqual(rtools.strand_confidence_threshold, 0.567)
        self.assertTrue(rtools._use_strand_correction)
