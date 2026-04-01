import unittest
from tempfile import NamedTemporaryFile
import os
from reditools.alignment_file import RTAlignmentFile
from reditools.region import Region
from .sam_gen import SAM, Sequence


class TestRTAlignmentFile(unittest.TestCase):
    def setUp(self):
        self.sam_obj = SAM()
        self.sam_obj.add_contig('chr1')
        self.refseq = self.sam_obj.genome['chr1']

        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8',
                                suffix='.fa') as f:
            self.genome_fname = f.name

        with NamedTemporaryFile(delete=False,
                                mode='w',
                                encoding='utf-8',
                                suffix='.bam') as f:
            self.bam_fname = f.name

    def tearDown(self):
        os.remove(self.genome_fname)
        os.remove(self.bam_fname)

    def test_fetch(self):
        self.sam_obj.add_read('chr1', Sequence(self.refseq[:40], 0))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[20:], 20))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[20:], 20))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[40:], 40))

        self.sam_obj.add_contig('chr2')
        self.sam_obj.add_read('chr2', Sequence(self.sam_obj.genome['chr2'], 0))

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        chr1_region = Region.from_string('chr1:60', self.bam_fname)
        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = list(rtaf.fetch())
            self.assertEqual(len(reads), 5)

            reads = list(rtaf.fetch(region='chr2'))
            self.assertEqual(len(reads), 1)

            reads = list(rtaf.fetch(region=chr1_region))
            self.assertEqual(len(reads), 3)

    def test_fetch_by_position(self):
        self.sam_obj.add_read('chr1', Sequence(self.refseq[:-20], 0))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[20:], 20))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[20:], 20))
        self.sam_obj.add_read('chr1', Sequence(self.refseq[40:], 40))

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

        with RTAlignmentFile(self.bam_fname) as rtaf:
            rtaf.exclude_reads = ['exclude_me']
            reads = list(rtaf.fetch())
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
            reads = list(rtaf.fetch())
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
            reads = list(rtaf.fetch())
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
            reads = list(rtaf.fetch())
            self.assertEqual(len(reads), 2)
            self.assertTrue(all(_.qname.startswith('se_good') for _ in reads))

    def test_check_pe_flags(self):
        for idx, flag in enumerate([83, 99]):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                qname=f'pe_good_{idx}',
            )
            self.sam_obj.add_read('chr1', read)
            self.sam_obj.add_read('chr1', read.make_pair())
        bad_flags = [73, 89, 137, 153, 329, 339, 345, 355, 393, 409]
        for idx, flag in enumerate(bad_flags):
            read = Sequence(
                self.refseq,
                0,
                flag=flag,
                qname=f'pe_bad_{idx}',
            )
            self.sam_obj.add_read('chr1', read)
            self.sam_obj.add_read('chr1', read.make_pair())

        self.sam_obj.genome.save_to_fasta(self.genome_fname)
        self.sam_obj.save_to_sam(self.bam_fname, self.genome_fname)

        with RTAlignmentFile(self.bam_fname) as rtaf:
            reads = list(rtaf.fetch())
            self.assertEqual(len(reads), 4)
            self.assertTrue(all(_.qname.startswith('pe_good') for _ in reads))
