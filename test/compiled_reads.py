import os
import unittest
from test.sam_gen import SAM, Sequence, ntf

from pysam import AlignmentFile

from reditools.compiled_reads import CompiledReads, RefFetch


class TestCompiledReads(unittest.TestCase):
    def setUp(self):
        self.fasta_fname = ntf(suffix='.fa')
        self.bam_fname = ntf(suffix='.bam')

    def tearDown(self):
        os.remove(self.fasta_fname)
        os.remove(self.bam_fname)

    def test_ref_seq_spliced(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=60)
        spliceseq = sam_obj.genome['chr1']
        spliceseq = spliceseq[:20] + spliceseq[40:60]
        sam_obj.add_read(
            'chr1',
            Sequence(spliceseq, 0, _cigar_str='20M20D20M'),
        )
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual(''.join(md_ref_fetch.get_refseq(read)), spliceseq)
        self.assertEqual(''.join(fa_ref_fetch.get_refseq(read)), spliceseq)

    def test_ref_seq_unspliced(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=60)
        refseq = sam_obj.genome['chr1']
        sam_obj.add_read('chr1', Sequence(refseq, 0))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual(''.join(md_ref_fetch.get_refseq(read)), refseq)
        self.assertEqual(''.join(fa_ref_fetch.get_refseq(read)), refseq)

    def test_ref_seq_snp(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=60)
        snpseq = list(sam_obj.genome['chr1'])
        snpseq[30] = 'A' if snpseq[30] == 'T' else 'T'
        snpseq = ''.join(snpseq)
        sam_obj.add_read('chr1', Sequence(snpseq, 0, _cigar_str='30M1X29M'))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        md_ref_fetch = RefFetch()
        fa_ref_fetch = RefFetch(self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
        self.assertEqual(
            ''.join(fa_ref_fetch.get_refseq(read)),
            sam_obj.genome['chr1'],
        )
        self.assertEqual(
            ''.join(md_ref_fetch.get_refseq(read)),
            sam_obj.genome['chr1'],
        )

    def test_se_strands(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1')
        ref_seq = sam_obj.genome['chr1']
        sam_obj.add_read(
            'chr1',
            Sequence(ref_seq, 0, flag=0, qname='read1'),
        )
        sam_obj.add_read(
            'chr1',
            Sequence(ref_seq[1:], 1, flag=16, qname='read2'),
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

    def test_pe_strands(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1')
        ref_seq = sam_obj.genome['chr1']
        sam_obj.add_read_pair('chr1', Sequence(ref_seq, 0, flag=99))
        sam_obj.add_read_pair('chr1', Sequence(ref_seq[1:], 1, flag=83))
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

    def test_trim(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=20)
        sam_obj.add_read('chr1', Sequence(sam_obj.genome['chr1'], 0))
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
            cr = CompiledReads(min_base_position=5, max_base_position=5)
            cr.add_reads([read])
            self.assertEqual(min(cr._nucleotides.keys()), 5)
            self.assertEqual(max(cr._nucleotides.keys()), 15)

    def test_base_quality(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=20)
        read = Sequence(sam_obj.genome['chr1'], 0, phred=range(20))
        sam_obj.add_read('chr1', read)
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read = next(af.fetch())
            cr = CompiledReads(min_base_quality=10)
            for _, _, phred, _ in cr._prep_read(read):
                self.assertTrue(phred >= 10)
            cr.add_reads([read])
            self.assertEqual(len(cr._nucleotides), 10)
