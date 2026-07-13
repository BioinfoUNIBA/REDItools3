from __future__ import annotations

import random
import unittest
from itertools import chain
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.fasta_file import MissingContigError, RTFastaFile


class TestRTFastaFile(unittest.TestCase):
    """Test cases for RTFastaFaile class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.contig1 = "test1"
        self.seq1 = self.random_seq(80)
        self.contig2 = "chrtest2"
        self.seq2 = self.random_seq(80)
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8",
        ) as stream:
            self.fasta_fname = stream.name
            stream.write(f">{self.contig1}\n{self.seq1}\n")
            stream.write(f">{self.contig2}\n{self.seq2}\n")

    def tearDown(self) -> None:
        """Post-check cleanup."""
        Path(self.fasta_fname).unlink()

    def test_get_base(self) -> None:
        """Check get_base() function."""
        with RTFastaFile(self.fasta_fname) as rff:
            positions = list(range(len(self.seq1)))
            fasta_seq = rff.get_base(self.contig1, *positions)
            self.assertEqual(self.seq1, "".join(fasta_seq))

    def test_get_base_splice(self) -> None:
        """Check get_base() function with spliced reads."""
        with RTFastaFile(self.fasta_fname) as rff:
            positions = list(chain(
                range(20),
                range(len(self.seq1) - 20, len(self.seq1)),
            ))
            fasta_seq = rff.get_base(self.contig1, *positions)
            self.assertEqual(
                self.seq1[:20] + self.seq1[-20:],
                "".join(fasta_seq),
            )

    def test_get_base_prefix(self) -> None:
        """Check get_base() function with contig name variants.

        Specifically, get_base should work regardless of whether a chromosome
        starts with "chr".
        """
        with RTFastaFile(self.fasta_fname) as rff:
            positions = range(len(self.seq2))
            fasta_seq = rff.get_base("test2", *positions)
            self.assertEqual(self.seq2, "".join(fasta_seq))

            fasta_seq = rff.get_base("chrtest2", *positions)
            self.assertEqual(self.seq2, "".join(fasta_seq))

    def test_get_base_missing_contig(self) -> None:
        """Check errors when accessing non-existent chromosomes."""
        with RTFastaFile(self.fasta_fname) as rff:
            with self.assertRaises(MissingContigError):
                rff.get_base("test3", 0)
            with self.assertRaises(MissingContigError):
                rff.get_base("chrtest3", 0)

    def test_get_base_out_of_bounds(self) -> None:
        """Check errors when accessing bases outside of chromosome boundns."""
        with RTFastaFile(self.fasta_fname) as rff, \
                self.assertRaises(IndexError):
            start = len(self.seq1) - 20
            stop = len(self.seq1) + 20
            positions = range(start, stop)
            seq_iter = rff.get_base(
                self.contig1,
                *positions,
            )
            list(seq_iter)

    @classmethod
    def random_seq(cls, length: int) -> str:
        """Generate a random nucleotide sequence.

        Parameters
        ----------
        length : int
            Sequence length.

        Returns
        -------
        str
            Random sequence.
        """
        sequence = [random.choice("ACTG") for _ in range(length)]
        return "".join(sequence)

