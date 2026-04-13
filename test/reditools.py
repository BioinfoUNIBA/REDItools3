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
        bam_file = ntf(suffix='.bam')
        fa_file = ntf(suffix='.fa')

        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=5)
        sam_obj.add_read('chr1', Sequence(sam_obj.genome['chr1'], 0))

        sam_obj.genome.save_to_fasta(fa_file)
        sam_obj.save_to_sam(bam_file, fa_file)

        rtam = AlignmentManager()
        rtam.add_file(bam_file)

        rtresult = next(
            self.rtools.analyze(rtam, Region.from_string('chr1', bam_file)),
        )
        self.assertEqual(rtresult.reference, sam_obj.genome['chr1'][0])

        new_genome = Genome()
        new_genome.add_contig(
            'chr1',
            sequence=''.join(
                self.complement[_] for _ in sam_obj.genome['chr1']
            ),
        )
        new_genome.save_to_fasta(fa_file)

        rtresult = next(
            self.rtools.analyze(rtam, Region.from_string('chr1', bam_file)),
        )
        self.assertEqual(rtresult.reference, sam_obj.genome['chr1'][0])

        os.remove(bam_file)
        os.remove(fa_file)

'''
    def __init__(self):
    def analyze(self, alignment_manager, region):
    def add_reference(self, reference_fname):
'''
