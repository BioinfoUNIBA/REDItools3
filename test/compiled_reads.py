"""Test cases for CompiledReads class."""
from __future__ import annotations

import unittest
from pathlib import Path
from test.sam_gen import SAM, Sequence, ntf

from pysam import AlignmentFile

from reditools.compiled_reads import CompiledReads, RefFetch


class TestCompiledReads(unittest.TestCase):
    """Test cases for CompiledReads class."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.fasta_fname = ntf(suffix=".fa")
        self.bam_fname = ntf(suffix=".bam")

    def tearDown(self) -> None:
        """Post-check cleanup."""
        Path(self.fasta_fname).unlink()
        Path(self.bam_fname).unlink()

    def test_ref_seq_spliced(self) -> None:
        """Check ability to get reference sequence from spliced reads."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=60)
        spliceseq = sam_obj.genome["chr1"]
        spliceseq = spliceseq[:20] + spliceseq[40:60]
        sam_obj.add_read(
            "chr1",
            Sequence(spliceseq, 0, _cigar_str="20M20D20M"),
        )
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual("".join(md_ref_fetch.get_refseq(read)), spliceseq)
        self.assertEqual("".join(fa_ref_fetch.get_refseq(read)), spliceseq)

    def test_ref_seq_unspliced(self) -> None:
        """Check ability to get reference sequence from contiguous reads."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=60)
        refseq = sam_obj.genome["chr1"]
        sam_obj.add_read("chr1", Sequence(refseq, 0))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual("".join(md_ref_fetch.get_refseq(read)), refseq)
        self.assertEqual("".join(fa_ref_fetch.get_refseq(read)), refseq)

    def test_ref_seq_snp(self) -> None:
        """Check ability to get reference sequence from reads with SNPs."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=60)
        snpseq_list = list(sam_obj.genome["chr1"])
        snpseq_list[30] = "A" if snpseq_list[30] == "T" else "T"
        sam_obj.add_read("chr1", Sequence(
            "".join(snpseq_list),
            0,
            _cigar_str="30M1X29M",
        ))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual(
            "".join(fa_ref_fetch.get_refseq(read)),
            sam_obj.genome["chr1"],
        )
        self.assertEqual(
            "".join(md_ref_fetch.get_refseq(read)),
            sam_obj.genome["chr1"],
        )

    def test_se_strands(self) -> None:
        """Check ability to recognize single-stranded data.

        Checks for unstranded, forward, and reverse strand data.
        """
        sam_obj = SAM()
        sam_obj.add_contig("chr1")
        ref_seq = sam_obj.genome["chr1"]
        sam_obj.add_read(
            "chr1",
            Sequence(ref_seq, 0, flag=0, read_name="read1"),
        )
        sam_obj.add_read(
            "chr1",
            Sequence(ref_seq[1:], 1, flag=16, read_name="read2"),
        )

        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            reads = list(af.fetch())

        cr = CompiledReads(strand=0)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [2, 2],
        )

        cr = CompiledReads(strand=1)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [True, False],
        )

        cr = CompiledReads(strand=2)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [False, True],
        )

    def test_pe_strands(self) -> None:
        """Check ability to recognize paired-end data.

        Checks for unstranded, forward, and reverse strand data.
        """
        sam_obj = SAM()
        sam_obj.add_contig("chr1")
        ref_seq = sam_obj.genome["chr1"]
        sam_obj.add_read_pair("chr1", Sequence(ref_seq, 0, flag=99))
        sam_obj.add_read_pair("chr1", Sequence(ref_seq[1:], 1, flag=83))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            reads = list(af.fetch())

        cr = CompiledReads(strand=0)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [2, 2, 2, 2],
        )

        cr = CompiledReads(strand=1)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [True, True, False, False],
        )

        cr = CompiledReads(strand=2)
        self.assertEqual(
            [cr.get_strand(_) for _ in reads],
            [False, False, True, True],
        )

    def test_trim(self) -> None:
        """Check read trimming."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=20)
        sam_obj.add_read("chr1", Sequence(sam_obj.genome["chr1"], 0))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
            cr = CompiledReads(min_base_position=5, max_base_position=5)
            cr.add_reads([read])
            self.assertEqual(min(cr._nucleotides.keys()), 5)
            self.assertEqual(max(cr._nucleotides.keys()), 15)

    def test_base_quality(self) -> None:
        """Check base quality filters."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=20)
        sam_obj.add_read("chr1", Sequence(
            sam_obj.genome["chr1"],
            0,
            phred_list=list(range(20)),
        ))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            algn_seg = next(af.fetch())
            mbq = 10
            cr = CompiledReads(min_base_quality=mbq)
            for _, _, phred, _ in cr._prep_read(algn_seg):
                self.assertTrue(phred >= mbq)
            cr.add_reads([algn_seg])
            self.assertEqual(len(cr._nucleotides), 10)

    def test_pop_range(self) -> None:
        """Check pop_range function."""
        sam_obj = SAM()
        sam_obj.add_contig("chr1", length=20)
        read = Sequence(
            sam_obj.genome["chr1"],
            0,
            phred_list=list(range(20)),
        )
        sam_obj.add_read("chr1", read)
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            cr = CompiledReads(min_base_quality=10)
            cr.add_reads([next(af.fetch())])
            cp_list = list(cr.pop_range(18, 30))
            self.assertEqual(len(cp_list), 2)

            self.assertEqual(cp_list[0].ref, sam_obj.genome["chr1"][18])
