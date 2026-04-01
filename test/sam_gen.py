import random
from Bio.Align import PairwiseAligner
import pysam.samtools as samtools
import os
from tempfile import NamedTemporaryFile
import re

aligner = PairwiseAligner(
    mismatch_score=-1,
    query_internal_open_gap_score=-1,
)


class Genome:
    def __init__(self):
        self.contigs = {}

    def _random_seq(self, length):
        return ''.join([random.choice('ACTG') for _ in range(length)])

    def __getitem__(self, contig_name):
        return self.contigs.get(contig_name, None)

    def add_contig(self, name=None, length=120, sequence=None):
        if name is None:
            n_contigs = len(self.contigs)
            name = f'contig{len(self.contigs)}'
            while name in self.contigs:
                n_contigs += 1
                name = f'contig{len(self.contigs)}'
        if sequence is None:
            self.contigs[name] = self._random_seq(length)
        else:
            self.contigs[name] = sequence
        return name

    def save_to_fasta(self, filename):
        with open(filename, 'w') as stream:
            for n, (name, sequence) in enumerate(self.contigs.items(), 1):
                stream.write(f'>{name} {n}\n{sequence}\n')
        samtools.faidx(filename)


class Sequence:
    read_n = 0

    def __init__(self, seq, start, flag=0, phred=None, mapq=None,
                 cigar_str=None, qname=None, pnext=0):
        self.seq = list(seq)
        self.start = start
        self.flag = flag
        self._cigar_str = cigar_str
        self.pnext = pnext
        if qname is None:
            self.qname = f'read{Sequence.read_n}'
            Sequence.read_n += 1
        else:
            self.qname = qname
        if phred is None:
            self.phred = [30 for _ in range(len(self.seq))]
        else:
            self.phred = phred
        self.mapq = mapq if mapq is not None else 255

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return ''.join(self.seq)

    def tlen(self, ref_seq):
        cigar = self.cigar_str(ref_seq)
        tlen = 0
        for count, op in re.findall(r'(?P<count>\d+)(?P<op>[A-Z])', cigar):
            count = int(count)
            if op not in ('S', 'I'):
                tlen += count
        if self.flag & 16:
            return -tlen
        return tlen

    def cigar_str(self, ref_seq):
        if self._cigar_str is not None:
            return self._cigar_str
        alignment = aligner.align(
            ref_seq[self.start:self.start + len(self)],
            str(self),
        )[0]
        cigar_pieces = []
        last_op = None
        op_n = 0
        for ref, query in zip(alignment[0], alignment[1]):
            if ref == '-':
                cigar_op = 'I'
            elif query == '-':
                if len(cigar_pieces) == 0:
                    cigar_op = 'S'
                else:
                    cigar_op = 'D'
            else:
                cigar_op = 'M'

            if cigar_op == last_op:
                op_n += 1
            else:
                if last_op is not None:
                    cigar_pieces.append(f'{op_n}{last_op}')
                last_op = cigar_op
                op_n = 1
        if last_op == 'D':
            last_op = 'S'
        cigar_pieces.append(f'{op_n}{last_op}')
        self._cigar_str = ''.join(cigar_pieces)
        return self._cigar_str

    def make_pair(self):
        return Sequence(
            seq=self.seq,
            start=self.start,
            flag=self.flag ^ 240,
            phred=self.phred,
            mapq=self.mapq,
            cigar_str=self._cigar_str,
            qname=self.qname,
            pnext=self.start,
        )
        self.pnext = self.start


class SAM:
    def __init__(self):
        self.genome = Genome()
        self.reads = {}

    def header(self):
        header = ['@HD\tVN:1.5']
        for contig_name, seq in self.genome.contigs.items():
            header.append(f'@SQ\tSN:{contig_name}\tLN:{len(seq)}')
        header.append(
            '@RG\tID:1\tSM:1_AAAAA\tLB:default\tPU:xxx.1\tPL:ILLUMINA',
        )
        header.append('@PG\tID:reditools\tPN:reditools\tCL:gen_sam.py')
        return '\n'.join(header)

    def add_contig(self, contig_name=None, length=120, sequence=None):
        contig_name = self.genome.add_contig(contig_name, length, sequence)
        self.reads[contig_name] = []
        return contig_name

    def add_read(self, contig_name, sequence_obj):
        self.reads[contig_name].append(sequence_obj)

    def __getitem__(self, contig_name):
        return self.reads[contig_name]

    def _covered_seqs(self, contig_name, position):
        return [
            idx for idx, seq in enumerate(self[contig_name])
            if seq.start <= position < seq.stop
        ]

    def sam_entries(self):
        n = 0
        for contig, reads in self.reads.items():
            ref_seq = self.genome[contig]
            for idx, sequence in enumerate(reads):
                yield "\t".join([str(_) for _ in [
                    sequence.qname,
                    sequence.flag,
                    contig,
                    sequence.start + 1,
                    sequence.mapq,
                    sequence.cigar_str(ref_seq),
                    ('*', '=')[sequence.flag & 1],
                    sequence.pnext + 1,
                    sequence.tlen(ref_seq),
                    str(sequence),
                    ''.join([chr(33 + _) for _ in sequence.phred]),
                ]])
                n += 1

    def save_to_sam(self, bam_filename, genome_filename):
        with NamedTemporaryFile(delete=False, mode='w', dir='.') as stream:
            sam_filename = stream.name
            stream.write(self.header())
            stream.write('\n')
            stream.write("\n".join(self.sam_entries()))
        md_sam = samtools.calmd(
            sam_filename,
            genome_filename,
            catch_stdout=True,
        )
        with open(sam_filename, 'w') as stream:
            stream.write(md_sam)
        samtools.sort('-o', bam_filename, sam_filename)
        samtools.index(bam_filename)
        os.remove(sam_filename)
