import unittest
from tempfile import NamedTemporaryFile
import os
from reditools.alignment_manager import AlignmentManager
from .sam_gen import SAM, Sequence


class TestAlignmentManager(unittest.TestCase):

    def test_propagation(self):
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.fa') as f:
            genome_fname = f.name

        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.bam') as f:
            bam_fname = f.name

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
        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.fa') as f:
            genome_fname = f.name

        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.bam') as f:
            bam_fname_1 = f.name

        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=80)
        refseq = sam_obj.genome['chr1']
        sam_obj.genome.save_to_fasta(genome_fname)

        sam_obj.add_read('chr1', Sequence(refseq, 0, qname='1_1'))
        sam_obj.add_read('chr1', Sequence(refseq[20:], 20, qname='1_2'))
        sam_obj.add_read('chr1', Sequence(refseq[40:], 40, qname='1_3'))
        sam_obj.save_to_sam(bam_fname_1, genome_fname)

        with NamedTemporaryFile(delete=False,
                                mode='w',
                                suffix='.bam') as f:
            bam_fname_2 = f.name

        sam_obj = SAM()
        sam_obj.add_contig('chr1', sequence=refseq)
        sam_obj.add_read('chr1', Sequence(refseq[20:], 20, qname='2_1'))
        sam_obj.add_read('chr1', Sequence(refseq[50:], 50, qname='2_2'))
        sam_obj.save_to_sam(bam_fname_2, genome_fname)

        rtam = AlignmentManager(min_length=10, min_quality=30)
        rtam.add_file(bam_fname_1)
        rtam.add_file(bam_fname_2)
        read_group_iter = rtam.fetch_by_position('chr1')

        first_group = next(read_group_iter)
        self.assertEqual(len(first_group), 1)
        self.assertEqual(first_group[0].qname, '1_1')

        second_group = next(read_group_iter)
        self.assertEqual(len(second_group), 2)
        self.assertIn('2_1', (_.qname for _ in second_group))
        self.assertIn('1_2', (_.qname for _ in second_group))

        os.remove(genome_fname)
        os.remove(bam_fname_1)
        os.remove(bam_fname_2)
