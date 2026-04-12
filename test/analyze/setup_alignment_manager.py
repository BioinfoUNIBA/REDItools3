import os
import unittest
from test.sam_gen import SAM, ntf

from reditools.region import Region
from reditools.tools.analyze.setup_alignment_manager import \
    setup_alignment_manager


class TestSetupAlignmentManager(unittest.TestCase):
    def test_setup(self):
        fasta_fname = ntf(suffix='.fa')
        bam_fname = ntf(suffix='.bam')

        sam_obj = SAM()
        sam_obj.add_contig('chr1', length=120)
        sam_obj.add_contig('chr2', length=80)
        sam_obj.add_contig('chr3', length=60)

        sam_obj.genome.save_to_fasta(fasta_fname)
        sam_obj.save_to_sam(bam_fname, fasta_fname)


        exclusions_fname = ntf(suffix='.bed')
        with open(exclusions_fname, 'w') as stream:
            stream.write('chr1\t100\t200\n')

        rtam = setup_alignment_manager(
            [bam_fname],
            min_read_quality=50,
            min_read_length=123,
            exclusions_file=exclusions_fname,
        )

        self.assertEqual(rtam.min_quality, 50)
        self.assertEqual(rtam.min_length, 123)
        self.assertIn(Region('chr1', 100, 200), rtam.exclude_set)

        os.remove(fasta_fname)
        os.remove(bam_fname)
        os.remove(exclusions_fname)
