import os
import unittest
from test.sam_gen import SAM, Genome, Sequence, ntf

from reditools.alignment_manager import AlignmentManager
from reditools.compiled_position import CompiledPosition
from reditools.reditools import REDItools
from reditools.region import Region


class TestREDItools(unittest.TestCase):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def setUp(self):
        self.rtools = REDItools()

        self.bam_file = ntf(suffix='.bam')
        self.fa_file = ntf(suffix='.fa')

        self.sam_obj = SAM()
        self.sam_obj.add_contig('chr1', length=10)
        self.sam_obj.add_read('chr1', Sequence(self.sam_obj.genome['chr1'], 0))

        self.sam_obj.genome.save_to_fasta(self.fa_file)
        self.sam_obj.save_to_sam(self.bam_file, self.fa_file)

        self.rtam = AlignmentManager()
        self.rtam.add_file(self.bam_file)

    def tearDown(self):
        os.remove(self.bam_file)
        os.remove(self.fa_file)

    def test_process_bases(self):
        cp = CompiledPosition(ref='A', position=1, contig='chr1')
        cp.add_base(30, '-', 'A')
        cp.add_base(30, '-', 'A')
        cp.add_base(30, '-', 'A')
        cp.add_base(30, '+', 'G')
        cp.add_base(30, '+', 'G')

        self.rtools.strand = 2
        rtresult = self.rtools._process_bases(cp)
        self.assertEqual(rtresult.reference, 'A')
        self.assertEqual(rtresult.strand, '*')
        self.assertEqual(rtresult.variants, ['AG'])

        self.rtools.strand = 1
        self.rtools.strand_confidence_threshold = 0.5
        rtresult = self.rtools._process_bases(cp)
        self.assertEqual(rtresult.strand, '-')
        self.assertEqual(rtresult.reference, 'A')
        self.assertEqual(rtresult.variants, ['AG'])

        self.rtools.use_strand_correction()
        rtresult = self.rtools._process_bases(cp)
        self.assertEqual(rtresult.strand, '-')
        self.assertEqual(rtresult.reference, 'T')
        self.assertEqual(rtresult.variants, [])

    def test_add_reference(self):
        rtresult = next(
            self.rtools.analyze(
                self.rtam,
                Region.from_string('chr1', self.bam_file),
            ),
        )
        self.assertEqual(rtresult.reference, self.sam_obj.genome['chr1'][0])

        new_genome = Genome()
        new_genome.add_contig(
            'chr1',
            sequence=''.join(
                self.complement[_] for _ in self.sam_obj.genome['chr1']
            ),
        )
        new_genome.save_to_fasta(self.fa_file)

        rtresult = next(
            self.rtools.analyze(
                self.rtam,
                Region.from_string('chr1', self.bam_file),
            ),
        )
        self.assertEqual(rtresult.reference, self.sam_obj.genome['chr1'][0])

    def test_region(self):
        rtresults = list(self.rtools.analyze(
            self.rtam,
            Region.from_string('chr1:3-7', self.bam_file),
        ))
        self.assertEqual(len(rtresults), 5)

        rtresults = list(self.rtools.analyze(
            self.rtam,
            Region.from_string('chr1:8-20', self.bam_file),
        ))
        self.assertEqual(len(rtresults), 3)
