import unittest
import os
from tempfile import NamedTemporaryFile
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.compiled_position import CompiledPosition
from reditools.region import Region

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __getattr__(self, attr):
        return self.__dict__.get(attr, None)

class TestRTChecks(unittest.TestCase):
    def setUp(self):
        self.bases = CompiledPosition(contig='chr1', position=1, ref='A')
        self.options = Namespace(
            max_editing_nucleotides=4,
            min_read_depth=0,
            min_edits=0,
            min_edits_per_nucleotide=0,
            variants=['all'],
        )

    def test_check_column_edit_frequency(self):
        self.options.min_edits = 1
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertIsNotNone(rtc.check(self.bases))

        self.options.min_edits = 0
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        self.options.min_edits = 1
        rtc = RTChecks(self.options)
        self.bases.add_base(quality=30, base='T', strand='*')
        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.options.min_edits = 3
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

    def test_check_column_min_edits(self):
        self.options.min_edits_per_nucleotide = 1
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='A', strand='*')
        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.options.min_edits_per_nucleotide = 2
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertIsNotNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertIsNotNone(rtc.check(self.bases))

    def test_check_min_read_depth(self):
        self.options.min_read_depth = 2
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertIsNotNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertIsNone(rtc.check(self.bases))

    def test_check_exclusions(self):
        with NamedTemporaryFile(delete=False, suffix='.bed', mode='w+') as stream:
            stream.write('chr1\t20\t30\n')
            bed_file = stream.name

        self.options.exclude_regions = [bed_file]

        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        with open(bed_file, mode='a') as stream:
            stream.write('chr1\t0\t10\n')
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

        with open(bed_file, mode='a') as stream:
            stream.write('chr2\t0\t10\n')
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

        os.remove(bed_file)

    def test_check_max_editing_nucleotides(self):
        self.options.max_editing_nucleotides = 1
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='A', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='T', strand='*')
        self.bases.add_base(quality=30, base='T', strand='*')
        self.assertIsNone(rtc.check(self.bases))

        self.bases.add_base(quality=30, base='C', strand='*')
        self.assertIsNotNone(rtc.check(self.bases))

        self.options.max_editing_nucleotides = 2
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

    def test_check_target_positions(self):
        with NamedTemporaryFile(delete=False, suffix='.bed', mode='w+') as stream:
            stream.write('chr1\t10\t20\n')
            bed_file = stream.name
        self.options.bed_file = [bed_file]
        rtc = RTChecks(self.options)
        self.assertIsNotNone(rtc.check(self.bases))

        with open(bed_file, mode='a') as stream:
            stream.write('chr1\t0\t20\n')
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        with open(bed_file, mode='a') as stream:
            stream.write('chr2\t0\t20\n')
        rtc = RTChecks(self.options)
        self.assertIsNone(rtc.check(self.bases))

        os.remove(bed_file)
