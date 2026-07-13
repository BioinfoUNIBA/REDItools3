"""Test cases for RTFastaFaile class."""
from __future__ import annotations

import random
import unittest
from itertools import chain
from pathlib import Path
from tempfile import NamedTemporaryFile
from test.sam_gen import Genome

from reditools.constants import bases
from reditools.fasta_file import (
    MissingContigError,
    PastContigEndError,
    RTFastaFile,
)


class TestRTFastaFile(unittest.TestCase):
    """Test cases for RTFastaFaile class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.genome = Genome()
        self.naked_contig_name = "test1"
        self.genome.add_contig(self.naked_contig_name, 80)
        self.chr_contig_name = "chrtest2"
        self.genome.add_contig(self.chr_contig_name, 80)
        with NamedTemporaryFile(
                delete=False,
                mode="w",
                encoding="utf-8",
        ) as stream:
            self.fasta_fname = stream.name
        self.genome.save_to_fasta(self.fasta_fname)

    def tearDown(self) -> None:
        """Post-check cleanup."""
        Path(self.fasta_fname).unlink()

    def test_get_base(self) -> None:
        """Check get_base() function."""
        refseq = self.genome[self.naked_contig_name]
        with RTFastaFile(self.fasta_fname) as rff:
            fasta_seq = rff.get_base(
                self.naked_contig_name,
                *range(len(refseq)),
            )
            self.assertEqual(refseq, "".join(fasta_seq))

    def test_get_base_splice(self) -> None:
        """Check get_base() function with spliced reads."""
        refseq = self.genome[self.naked_contig_name]
        with RTFastaFile(self.fasta_fname) as rff:
            positions = chain(
                range(20),
                range(len(refseq) - 20, len(refseq)),
            )
            fasta_seq = rff.get_base(self.naked_contig_name, *positions)
            self.assertEqual(
                refseq[:20] + refseq[-20:],
                "".join(fasta_seq),
            )

    def test_get_base_prefix(self) -> None:
        """Check get_base() function with contig name variants.

        Specifically, get_base should work regardless of whether a chromosome
        starts with "chr".
        """
        refseq = self.genome[self.chr_contig_name]
        with RTFastaFile(self.fasta_fname) as rff:
            fasta_seq = rff.get_base(
                self.chr_contig_name,
                *range(len(refseq)),
            )
            self.assertEqual(refseq, "".join(fasta_seq))

            fasta_seq = rff.get_base(
                f"chr{self.chr_contig_name}",
                *range(len(refseq)),
            )
            self.assertEqual(refseq, "".join(fasta_seq))

    def test_get_base_missing_contig(self) -> None:
        """Check errors when accessing non-existent chromosomes."""
        with RTFastaFile(self.fasta_fname) as rff:
            with self.assertRaises(MissingContigError):
                list(rff.get_base("test3", 0))
            with self.assertRaises(MissingContigError):
                list(rff.get_base("chrtest3", 0))

    def test_get_base_out_of_bounds(self) -> None:
        """Check errors when accessing bases outside of chromosome boundns."""
        refseq = self.genome[self.naked_contig_name]
        with RTFastaFile(self.fasta_fname) as rff, \
                self.assertRaises(PastContigEndError):
            positions = range(
                len(refseq) - 20,
                len(refseq) + 20,
            )
            seq_iter = rff.get_base(
                self.naked_contig_name,
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
        sequence = [random.choice(bases) for _ in range(length)]
        return "".join(sequence)
