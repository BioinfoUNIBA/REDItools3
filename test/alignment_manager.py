import unittest
import os
from reditools.alignment_manager import AlignmentManager
from test.sam_gen import SAM, Sequence, ntf


class TestAlignmentManager(unittest.TestCase):

    def test_propagation(self):
        genome_fname = ntf(suffix='.fa')
        bam_fname = ntf(suffix='.bam')

        sam_obj = SAM()
        sam_obj.add_contig('chr1')
        sam_obj.genome.save_to_fasta(genome_fname)
        sam_obj.save_to_sam(bam_fname, genome_fname)

        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(bam_fname)

        self.assertEqual(rtam._bams[0]._min_length, 10)
        self.assertEqual(rtam._bams[0]._min_quality, 30)

        os.remove(genome_fname)
        os.remove(bam_fname)

    def test_fetch_by_position(self):
        genome_fname, bam_fnames = self.setup_dummy_data()

        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(bam_fnames[0])
        rtam.add_file(bam_fnames[1])
        read_groups = list(rtam.fetch_by_position('chr1'))

        self.assertEqual(len(read_groups[0]), 1)
        self.assertEqual(read_groups[0][0].qname, '1_1')

        self.assertEqual(len(read_groups[1]), 2)
        self.assertIn('2_1', (_.qname for _ in read_groups[1]))
        self.assertIn('1_2', (_.qname for _ in read_groups[1]))

        os.remove(genome_fname)
        for fname in bam_fnames:
            os.remove(fname)

    def setup_dummy_data(self):
        genome_fname = ntf(suffix='.fa')
        bam_fnames = [ntf(suffix='.bam') for _ in range(2)]

        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=80)
        refseq = sam_obj.genome['chr1']
        sam_obj.genome.save_to_fasta(genome_fname)

        sam_obj.add_read('chr1', Sequence(refseq, 0, qname='1_1'))
        sam_obj.add_read('chr1', Sequence(refseq[20:], 20, qname='1_2'))
        sam_obj.add_read('chr1', Sequence(refseq[40:], 40, qname='1_3'))
        sam_obj.save_to_sam(bam_fnames[0], genome_fname)

        sam_obj = SAM()
        sam_obj.add_contig('chr1', sequence=refseq)
        sam_obj.add_read('chr1', Sequence(refseq[20:], 20, qname='2_1'))
        sam_obj.add_read('chr1', Sequence(refseq[50:], 50, qname='2_2'))
        sam_obj.save_to_sam(bam_fnames[1], genome_fname)

        return genome_fname, bam_fnames
