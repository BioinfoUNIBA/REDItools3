import unittest
import os
from reditools.fasta_file import RTFastaFile
from reditools.compiled_reads import CompiledReads
from pysam import AlignmentFile
from tempfile import NamedTemporaryFile
from .sam_gen import SAM, Sequence


class TestCompiledReads(unittest.TestCase):
    def setUp(self):
        with NamedTemporaryFile(delete=False, suffix='.fa') as f:
            self.fasta_fname = f.name
        with NamedTemporaryFile(delete=False, suffix='.bam') as f:
            self.bam_fname = f.name

    def tearDown(self):
        os.remove(self.fasta_fname)
        os.remove(self.bam_fname)

    def test_ref_seq(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=60)
        ref_seq = sam_obj.genome['chr1']
        normal_read = Sequence(ref_seq, 0, cigar_str='60M')
        splice_seq = ref_seq[0:20] + ref_seq[40:60]
        spliced_read = Sequence(splice_seq, 0, cigar_str='20M20D20M')
        sam_obj.add_read('chr1', normal_read)
        sam_obj.add_read('chr1', spliced_read)

        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        cr = CompiledReads()
        cr.add_reference(RTFastaFile(self.fasta_fname))

        with AlignmentFile(self.bam_fname) as af:
            for read, seq in zip(af.fetch(), [ref_seq, splice_seq]):
                self.assertEqual(
                    ''.join(cr._get_ref_from_read(read)),
                    seq,
                )
                self.assertEqual(
                    ''.join(cr._get_ref_from_fasta(read)),
                    seq,
                )

    def test_se_strands(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1')
        ref_seq = sam_obj.genome['chr1']
        read1 = Sequence(ref_seq, 0, flag=0, qname='read1')
        read2 = Sequence(ref_seq[1:], 1, flag=16, qname='read2')
        sam_obj.add_read('chr1', read1)
        sam_obj.add_read('chr1', read2)

        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            read1, read2 = af.fetch()

            cr = CompiledReads(strand=0)
            self.assertEqual(cr.get_strand(read1), 2)
            self.assertEqual(cr.get_strand(read2), 2)

            cr = CompiledReads(strand=1)
            self.assertEqual(cr.get_strand(read1), True)
            self.assertEqual(cr.get_strand(read2), False)

            cr = CompiledReads(strand=2)
            self.assertEqual(cr.get_strand(read1), False)
            self.assertEqual(cr.get_strand(read2), True)

    def test_pe_strands(self):
        sam_obj = SAM()
        sam_obj.add_contig('chr1')
        ref_seq = sam_obj.genome['chr1']
        read_83 = Sequence(ref_seq[1:], 1, flag=83)
        read_99 = Sequence(ref_seq, 0, flag=99)
        sam_obj.add_read('chr1', read_83)
        sam_obj.add_read('chr1', read_83.make_pair())
        sam_obj.add_read('chr1', read_99)
        sam_obj.add_read('chr1', read_99.make_pair())
        sam_obj.genome.save_to_fasta(self.fasta_fname)
        sam_obj.save_to_sam(self.bam_fname, self.fasta_fname)

        with AlignmentFile(self.bam_fname) as af:
            cr = CompiledReads(strand=0)
            for read in af.fetch():
                self.assertEqual(cr.get_strand(read), 2)

            cr = CompiledReads(strand=1)
            for read, strand in zip(af.fetch(), [True, True, False, False]):
                self.assertEqual(cr.get_strand(read), strand)

            cr = CompiledReads(strand=2)
            for read, strand in zip(af.fetch(), [False, False, True, True]):
                self.assertEqual(cr.get_strand(read), strand)

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
            for ref_pos, read_base, phred, ref_base in cr._prep_read(read):
                self.assertTrue(phred >= 10)
            cr.add_reads([read])
            self.assertEqual(len(cr._nucleotides), 10)
