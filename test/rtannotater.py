"""Test cases for RTAnnotater class."""
from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.compiled_position import RTResult
from reditools.constants import analyze_output_fieldnames
from reditools.rtannotater import AnalyzeMismatchError, RTAnnotater


class TestRTAnnotater(unittest.TestCase):
    """Test cases for RTAnnotater class."""

    def setUp(self) -> None:
        """Preflight setup."""
        self.rna = RTResult(
            reference="A",
            strand="*",
            position=1,
            contig="chr1",
            base_counts={"A": 1, "C": 0, "G": 3, "T": 0},
            variants=["AG"],
            depth=4,
            edit_ratio=0.25,
            mean_quality=30,
        )

        self.dna = RTResult(
            reference="A",
            strand="*",
            position=1,
            contig="chr1",
            base_counts={"A": 1, "C": 0, "G": 0, "T": 0},
            variants=[],
            depth=1,
            edit_ratio=0,
            mean_quality=35,
        )

    def test_cmp_position(self) -> None:
        """Check cmp_position() method."""
        contig_order = {
            "chrZ": 1,
            "chr1": 2,
            "chr2": 3,
        }

        self.rna.contig = "chr1"
        self.rna.position = 5

        self.dna.contig = "chr1"
        self.dna.position = 5

        rta = RTAnnotater(contig_order)
        self.assertEqual(rta.cmp_position(self.rna, self.dna), 0)

        self.rna.contig = "chrZ"
        self.dna.contig = "chr2"
        self.assertTrue(rta.cmp_position(self.rna, self.dna) < 0)

        self.rna.contig = "chr1"
        self.rna.position = 10

        self.dna.contig = "chr1"
        self.dna.position = 1
        self.assertTrue(rta.cmp_position(self.rna, self.dna) > 0)

    def test_annotate_complement_row(self) -> None:
        """Check annotate_row() on minus strand."""
        rta = RTAnnotater({}, True)  # noqa: FBT003

        self.rna.strand = "+"
        self.dna.strand = "-"
        self.dna.reference = "T"

        rta.annotate_row(self.rna, self.dna)
        output = self.rna.output_list()

        self.assertEqual(output[11], "[0, 0, 0, 1]")
        self.assertEqual(output[12], "-")

        self.dna.reference = "T"
        self.dna.variants = ["C"]
        rta.annotate_row(self.rna, self.dna)
        output = self.rna.output_list()
        self.assertEqual(output[12], "AG")

    def test_annotate_no_complement_row(self) -> None:
        """Check annotate_row() on minus strand without complementing."""
        rta = RTAnnotater({})

        self.rna.strand = "+"
        self.dna.strand = "-"
        self.dna.reference = "T"

        rta.annotate_row(self.rna, self.dna)
        output = self.rna.output_list()

        self.assertEqual(output[11], "[1, 0, 0, 0]")
        self.assertEqual(output[12], "-")

        self.dna.variants = ["C"]
        rta.annotate_row(self.rna, self.dna)
        output = self.rna.output_list()
        self.assertEqual(output[12], "TC")

    def test_mismatched_reference(self) -> None:
        """Check error handling for reference mismatch."""
        rta = RTAnnotater({})
        self.rna.reference = "A"
        self.dna.reference = "G"
        with self.assertRaises(AnalyzeMismatchError):
            rta.annotate_row(self.rna, self.dna)

    def test_merge_files(self) -> None:
        """Check merge_files() method."""
        with NamedTemporaryFile(
                delete=False,
                suffix=".out",
                mode="w",
        ) as stream:
            rna_file = stream.name
            stream.write("\t".join(analyze_output_fieldnames))
            stream.write("\n")
            stream.write("\t".join([
                "chr1", "3", "A", "*", "5", "30.00", "[5, 0, 0, 0]", "-",
                "0.00", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")
            stream.write("\t".join([
                "chr1", "5", "A", "*", "5", "30.00", "[2, 0, 3, 0]", "AG",
                "0.60", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")
            stream.write("\t".join([
                "chr2", "7", "G", "*", "5", "30.00", "[2, 0, 3, 0]", "GA",
                "0.40", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")

        with NamedTemporaryFile(
                delete=False,
                suffix=".out",
                mode="w",
        ) as stream:
            dna_file = stream.name
            stream.write("\t".join(analyze_output_fieldnames))
            stream.write("\n")
            stream.write("\t".join([
                "chrZ", "3", "A", "*", "5", "30.00", "[5, 0, 0, 0]", "-",
                "0.00", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")
            stream.write("\t".join([
                "chr1", "5", "A", "*", "5", "35.00", "[2, 0, 0, 3]", "AT",
                "0.60", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")
            stream.write("\t".join([
                "chr1", "9", "G", "*", "5", "30.00", "[2, 0, 3, 0]", "GA",
                "0.40", "-", "-", "-", "-", "-",
            ]))
            stream.write("\n")

        rta = RTAnnotater({"chrZ": 1, "chr1": 2, "chr2": 3})
        annotated_data = list(rta.merge_files(rna_file, dna_file))

        row = annotated_data.pop(0)
        self.assertEqual(row.output_list(), [
            "chr1", "3", "A", "*", "5", "30.00", "[5, 0, 0, 0]", "-", "0.00",
            "-", "-", "-", "-", "-",
        ])

        row = annotated_data.pop(0)
        self.assertEqual(row.output_list(), [
            "chr1", "5", "A", "*", "5", "30.00", "[2, 0, 3, 0]", "AG", "0.60",
            "5", "35.00", "[2, 0, 0, 3]", "AT", "0.60",
        ])

        row = annotated_data.pop(0)
        self.assertEqual(row.output_list(), [
            "chr2", "7", "G", "*", "5", "30.00", "[2, 0, 3, 0]", "GA", "0.40",
            "-", "-", "-", "-", "-",
        ])
        self.assertEqual(len(annotated_data), 0)
        Path(rna_file).unlink()
        Path(dna_file).unlink()
