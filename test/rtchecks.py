import unittest
from reditools import rtchecks
from reditools.compiled_position import CompiledPosition
from reditools.reditools import REDItools
from reditools.region import Region


class TestRTChecks(unittest.TestCase):
    def setUp(self):
        self.bases = CompiledPosition(contig='chr1', position=1, ref='A')
        self.rtools = REDItools()

    def test_check_column_edit_frequency(self):
        self.rtools.min_edits = 1
        self.assertFalse(
            rtchecks.check_column_edit_frequency(self.rtools, self.bases),
        )
        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertFalse(
            rtchecks.check_column_edit_frequency(self.rtools, self.bases),
        )
        self.rtools.min_edits = 0
        self.assertTrue(
            rtchecks.check_column_edit_frequency(self.rtools, self.bases),
        )

        self.rtools.min_edits = 1
        self.bases.add_base(quality=30, base='T', strand='*')
        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertTrue(
            rtchecks.check_column_edit_frequency(self.rtools, self.bases),
        )
        self.rtools.min_edits = 3
        self.assertFalse(
            rtchecks.check_column_edit_frequency(self.rtools, self.bases),
        )

    def test_check_column_min_edits(self):
        self.rtools.min_edits_per_nucleotide = 1
        self.assertTrue(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='A', strand='*')
        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertTrue(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

        self.rtools.min_edits_per_nucleotide = 2
        self.assertTrue(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertFalse(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertTrue(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertFalse(
            rtchecks.check_column_min_edits(self.rtools, self.bases),
        )

    def test_check_column_min_length(self):
        self.rtools.min_column_length = 1
        self.assertFalse(
            rtchecks.check_column_min_length(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertTrue(
            rtchecks.check_column_min_length(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertTrue(
            rtchecks.check_column_min_length(self.rtools, self.bases),
        )

    def test_check_column_quality(self):
        self.rtools.min_read_quality = 30
        self.assertFalse(
            rtchecks.check_column_quality(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertTrue(
            rtchecks.check_column_quality(self.rtools, self.bases),
        )

        self.bases.add_base(quality=14, base='C', strand='*')
        self.assertFalse(
            rtchecks.check_column_quality(self.rtools, self.bases),
        )

    def test_check_exclusions(self):
        self.assertTrue(
            rtchecks.check_exclusions(self.rtools, self.bases),
        )

        region = Region(contig='chr1', start=20, stop=30)
        self.rtools.add_exclude_regions([region])
        self.assertTrue(
            rtchecks.check_exclusions(self.rtools, self.bases),
        )

        region = Region(contig='chr1', start=0, stop=10)
        self.rtools.add_exclude_regions([region])
        self.assertFalse(
            rtchecks.check_exclusions(self.rtools, self.bases),
        )

        region = Region(contig='chr2', start=0, stop=10)
        self.rtools.add_exclude_regions([region])
        self.assertFalse(
            rtchecks.check_exclusions(self.rtools, self.bases),
        )

    def test_check_max_alts(self):
        self.rtools.max_alts = 1
        self.assertTrue(
            rtchecks.check_max_alts(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertTrue(
            rtchecks.check_max_alts(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='T', strand='*')
        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertTrue(
            rtchecks.check_max_alts(self.rtools, self.bases),
        )

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertFalse(
            rtchecks.check_max_alts(self.rtools, self.bases),
        )

        self.rtools.max_alts = 2
        self.assertTrue(
            rtchecks.check_max_alts(self.rtools, self.bases),
        )

    def test_check_target_positions(self):
        self.assertFalse(
            rtchecks.check_target_positions(self.rtools, self.bases),
        )

        region = Region(contig='chr1', start=10, stop=20)
        self.rtools.add_target_regions([region])
        self.assertFalse(
            rtchecks.check_target_positions(self.rtools, self.bases),
        )

        region = Region(contig='chr1', start=0, stop=20)
        self.rtools.add_target_regions([region])
        self.assertTrue(
            rtchecks.check_target_positions(self.rtools, self.bases),
        )

        region = Region(contig='chr2', start=0, stop=20)
        self.rtools.add_target_regions([region])
        self.assertTrue(
            rtchecks.check_target_positions(self.rtools, self.bases),
        )

    def test_rtchecks(self):
        rtc = rtchecks.RTChecks()
        true_fn = lambda a, b: True  # noqa: E731
        false_fn = lambda a, b: False  # noqa: E731

        self.assertTrue(rtc.check(self.rtools, self.bases))

        rtc.add(true_fn)
        self.assertTrue(rtc.check(self.rtools, self.bases))

        rtc.add(false_fn)
        self.assertFalse(rtc.check(self.rtools, self.bases))

        rtc.discard(false_fn)
        self.assertTrue(rtc.check(self.rtools, self.bases))
