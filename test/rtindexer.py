"""Test cases for RTIndexer class."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.rtindexer import RTIndexer


class TestRTIndexer(unittest.TestCase):
    """Test cases for RTIndexer class."""

    test_data = (
        {
            "Region": "chr1",
            "Position": 1,
            "Reference": "A",
            "BaseCount[A,C,G,T]": "[10, 0, 0, 0]",
        },
        {
            "Region": "chr1",
            "Position": 2,
            "Reference": "A",
            "BaseCount[A,C,G,T]": "[0, 0, 10, 0]",
        },
        {
            "Region": "chr1",
            "Position": 3,
            "Reference": "G",
            "BaseCount[A,C,G,T]": "[0, 10, 10, 0]",
        },
    )

    def setUp(self) -> None:
        """Pre-flight setup."""
        with NamedTemporaryFile(
                delete=False,
                suffix=".out",
                mode="wt",
        ) as stream:
            self.output_filename = stream.name
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    "Region",
                    "Position",
                    "Reference",
                    "BaseCount[A,C,G,T]",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(self.test_data)

        with NamedTemporaryFile(
                delete=False,
                suffix=".bed",
                mode="wt",
        ) as stream:
            self.bed_filename = stream.name
            stream.write("chr1\t0\t2\n")

    def tearDown(self) -> None:
        """Post-checks cleanup."""
        Path(self.output_filename).unlink()
        Path(self.bed_filename).unlink()

    def test_baseline(self) -> None:
        """Check calc_index() method."""
        rti = RTIndexer()
        rti.add_rt_output(self.output_filename)
        self.assertEqual(rti.calc_index(), {
            "A-C": 0,
            "A-T": 0,
            "A-G": 50,
            "C-A": 0,
            "C-T": 0,
            "C-G": 0,
            "G-A": 0,
            "G-C": 50,
            "G-T": 0,
            "T-A": 0,
            "T-C": 0,
            "T-G": 0,
        })

    def test_region(self) -> None:
        """Check do_ignore() method."""
        rti = RTIndexer(region=("chr1", 100, 200))
        self.assertFalse(rti.do_ignore({"Region": "chr1", "Position": "150"}))
        self.assertTrue(rti.do_ignore({"Region": "chr1", "Position": "50"}))
        self.assertTrue(rti.do_ignore({"Region": "chr1", "Position": "250"}))
        self.assertTrue(rti.do_ignore({"Region": "chr2", "Position": "150"}))

        rti =RTIndexer(region=("chr1", 100, None))
        self.assertFalse(rti.do_ignore({"Region": "chr1", "Position": "150"}))
        self.assertTrue(rti.do_ignore({"Region": "chr1", "Position": "50"}))
        self.assertTrue(rti.do_ignore({"Region": "chr2", "Position": "150"}))

    def test_targets(self) -> None:
        """Check add_target_from_bed() method."""
        rti = RTIndexer()
        rti.add_target_from_bed(self.bed_filename)
        self.assertFalse(rti.do_ignore({"Region": "chr1", "Position": "1"}))
        self.assertTrue(rti.do_ignore({"Region": "chr1", "Position": "2"}))
        self.assertTrue(rti.do_ignore({"Region": "chr2", "Position": "1"}))

    def test_exclusions(self) -> None:
        """Check add_exclusions_from_bed method."""
        rti = RTIndexer()
        rti.add_exclusions_from_bed(self.bed_filename)
        self.assertTrue(rti.do_ignore({"Region": "chr1", "Position": "1"}))
        self.assertFalse(rti.do_ignore({"Region": "chr1", "Position": "2"}))
        self.assertFalse(rti.do_ignore({"Region": "chr2", "Position": "1"}))
