import os
import unittest
from test.sam_gen import SAM, Sequence, ntf

from reditools.alignment_file import RTAlignmentFile


class TestRTAlignmentFile(unittest.TestCase):
    def setUp(self):
        self.sam_obj = SAM()
        self.sam_obj.add_contig('chr1', length=60)
        self.refseq = self.sam_obj.genome['chr1']

        self.genome_fname = ntf(suffix='.fa')
        self.bam_fname = ntf(suffix='.bam')

    def tearDown(self):
        os.remove(self.genome_fname)
        os.remove(self.bam_fname)

    def test_fetch_by_position(self):
        for start, stop in (
                (0, 20),
                (20, None),
                (20, None),
                (40, None),
        ):
            read_seq = self.refseq[start:stop]
            self.sam_obj.add_read('chr1', Sequence(read_seq, start))

        self.sam_obj.add_contig('chr2')
        self.sam_obj.add_read('chr2', Sequence(self.sam_obj.genome['chr2'], 0))

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads_iter = rtaf.fetch_by_position('chr1')
            self.assertEqual(len(next(reads_iter)), 1)
            self.assertEqual(len(next(reads_iter)), 2)
            self.assertEqual(len(next(reads_iter)), 1)

    def test_exclude_reads(self):
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq, 0, qname='exclude_me'),
        )
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq, 0, qname='include_me'),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(
                self.bam_fname,
                excluded_read_names={'exclude_me'},
        ) as rtaf:
            reads = next(rtaf.fetch_by_position('chr1'))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].qname, 'include_me')

    def test_check_quality(self):
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq, 0, mapq=10),
        )
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq, 0, mapq=30, qname='include_me'),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname, min_quality=20) as rtaf:
            reads = next(rtaf.fetch_by_position('chr1'))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].qname, 'include_me')

    def test_check_length(self):
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq[:20], 0),
        )
        self.sam_obj.add_read(
            'chr1',
            Sequence(self.refseq, 0, qname='include_me'),
        )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname, min_length=30) as rtaf:
            reads = next(rtaf.fetch_by_position('chr1'))
            self.assertEqual(len(reads), 1)
            self.assertEqual(reads[0].qname, 'include_me')

    def test_check_se_flags(self):
        for idx, flag in enumerate([0, 16]):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                qname=f'se_good_{idx}',
            )
            self.sam_obj.add_read('chr1', read)
        for idx, flag in enumerate([4, 256, 272, 512, 1024, 2048, 2064]):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                qname=f'se_bad_{idx}',
            )
            self.sam_obj.add_read('chr1', read)

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = next(rtaf.fetch_by_position('chr1'))
            self.assertEqual(len(reads), 2)
            self.assertTrue(all(_.qname.startswith('se_good') for _ in reads))

    def test_check_pe_flags(self):
        for idx, flag in enumerate([83, 99]):
            self.sam_obj.add_read_pair(
                'chr1',
                Sequence(
                    self.refseq,
                    0,
                    flag=flag,
                    qname=f'pe_good_{idx}',
                ),
            )
        bad_flags = [73, 89, 137, 153, 329, 339, 345, 355, 393, 409]
        for idx, flag in enumerate(bad_flags):
            self.sam_obj.add_read_pair(
                'chr1',
                Sequence(
                    self.refseq,
                    0,
                    flag=flag,
                    qname=f'pe_bad_{idx}',
                ),
            )

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = next(rtaf.fetch_by_position('chr1'))
            self.assertEqual(len(reads), 4)
            self.assertTrue(all(_.qname.startswith('pe_good') for _ in reads))
