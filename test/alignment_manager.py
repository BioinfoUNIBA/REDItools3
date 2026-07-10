from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, Sequence, ntf

from reditools.alignment_manager import AlignmentManager


class TestAlignmentManager(unittest.TestCase):

    def test_propagation(self) -> None:
        genome_fname = ntf(suffix=".fa")
        bam_fname = ntf(suffix=".bam")

        sam_obj = SAM()
        sam_obj.add_contig("chr1")
        sam_obj.genome.save_to_fasta(genome_fname)
        sam_obj.save_to_sam(bam_fname, genome_fname)

        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(bam_fname)

        self.assertEqual(rtam._bams[0].readqc.min_length, 10)
        self.assertEqual(rtam._bams[0].readqc.min_quality, 30)

        Path(genome_fname).unlink()
        Path(bam_fname).unlink()

    def test_fetch_by_position(self) -> None:
        genome_fname, bam_fnames = self.setup_dummy_data()

        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(bam_fnames[0])
        rtam.add_file(bam_fnames[1])

        read_iter = rtam.fetch_by_position("chr1")

        read_group = next(read_iter)
        self.assertEqual(len(read_group), 1)
        self.assertEqual(read_group[0].qname, "1_1")
        self.assertEqual(rtam.next_read_start, 20)

        read_group = next(read_iter)
        self.assertEqual(len(read_group), 2)
        self.assertIn("2_1", (_.qname for _ in read_group))
        self.assertIn("1_2", (_.qname for _ in read_group))
        self.assertEqual(rtam.next_read_start, 40)

        Path(genome_fname).unlink()
        for fname in bam_fnames:
            Path(fname).unlink()

    def setup_dummy_data(self) -> tuple[str, list[str]]:
        genome_fname = ntf(suffix=".fa")
        bam_fnames = [ntf(suffix=".bam") for _ in range(2)]

        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=80)
        refseq = sam_obj.genome["chr1"]
        sam_obj.genome.save_to_fasta(genome_fname)

        sam_obj.add_read("chr1", Sequence(refseq, 0, read_name="1_1"))
        sam_obj.add_read("chr1", Sequence(refseq[20:], 20, read_name="1_2"))
        sam_obj.add_read("chr1", Sequence(refseq[40:], 40, read_name="1_3"))
        sam_obj.save_to_sam(bam_fnames[0], genome_fname)

        sam_obj = SAM()
        sam_obj.add_contig("chr1", sequence=refseq)
        sam_obj.add_read("chr1", Sequence(refseq[20:], 20, read_name="2_1"))
        sam_obj.add_read("chr1", Sequence(refseq[50:], 50, read_name="2_2"))
        sam_obj.save_to_sam(bam_fnames[1], genome_fname)

        return genome_fname, bam_fnames
