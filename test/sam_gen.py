import random
from Bio.Align import PairwiseAligner
from pysam import samtools
import os
from tempfile import NamedTemporaryFile
import re
from dataclasses import dataclass, InitVar


class Genome:
    def __init__(self):
        self.contigs = {}

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
            for idx, (name, sequence) in enumerate(self.contigs.items(), 1):
                stream.write(f'>{name} {idx}\n{sequence}\n')
        samtools.faidx(filename)

    @classmethod
    def _random_seq(cls, length):
        return ''.join([random.choice('ACTG') for _ in range(length)])


@dataclass
class Sequence: 
    seq: str
    start: int
    flag: int = 0
    phred: InitVar[list | None] = None
    mapq: int = 255
    _cigar_str: str | None = None
    qname: InitVar[str | None] = None
    pnext: int = 0
    
    read_n = 0
    flag_reverse_strand = 16
    phred_default = 30
    aligner = PairwiseAligner(
        mismatch_score=-1,
        query_internal_open_gap_score=-1,
    )

    def __post_init__(self, phred, qname):
        if phred is None:
            self.phred = [self.phred_default for _ in range(len(self.seq))]
        else:
            self.phred = phred

        if qname is None:
            self.qname = self.next_read_name()
        else:
            self.qname = qname

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def tlen(self, ref_seq):
        cigar = self.cigar_str(ref_seq)
        tlen = 0
        for count, op in re.findall(r'(?P<count>\d+)(?P<op>[A-Z])', cigar):
            count = int(count)
            if op not in ('S', 'I'):
                tlen += count
        if self.flag & Sequence.flag_reverse_strand:
            return -tlen
        return tlen

    def cigar_str(self, ref_seq):
        if self._cigar_str is not None:
            return self._cigar_str
        alignment = Sequence.aligner.align(
            ref_seq[self.start:self.start + len(self)],
            str(self),
        )[0]
        cigar_iter = self.assemble_cigar_list(alignment[0], alignment[1])
        cigar_pieces = [f'{length}{op}' for length, op in cigar_iter]
        self._cigar_str = ''.join(cigar_pieces)
        return self._cigar_str

    def make_pair(self):
        return Sequence(
            seq=self.seq,
            start=self.start,
            flag=self.pair_flag(self.flag),
            phred=self.phred,
            mapq=self.mapq,
            _cigar_str=self._cigar_str,
            qname=self.qname,
            pnext=self.start,
        )

    @classmethod
    def pair_flag(cls, flag_value):
        return flag_value ^ 240

    @classmethod
    def cigar_op(cls, ref_base, query_base):
        if ref_base == '-':
            return 'I'
        if query_base == '-':
            return 'D'
        if query_base == ref_base:
            return 'M'
        return 'X'

    @classmethod
    def assemble_cigar_list(cls, algn_ref, algn_query):
        last_op = None
        op_n = 0
        for ref_base, query_base in zip(algn_ref, algn_query):
            cigar_op = cls.cigar_op(ref_base, query_base)
            if cigar_op == last_op:
                op_n += 1
            else:
                if last_op is not None:
                    yield (op_n, cigar_op)
                last_op = cigar_op
                op_n = 1
        yield (op_n, cigar_op)

    @classmethod
    def next_read_name(cls):
        cls.read_n += 1
        return f'read{cls.read_n}'


class SAM:
    def __init__(self):
        self.genome = Genome()
        self.reads = {}

    def __getitem__(self, contig_name):
        return self.reads[contig_name]

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

    def add_read_pair(self, contig_name, sequence_obj):
        self.add_read(contig_name, sequence_obj)
        self.add_read(contig_name, sequence_obj.make_pair())

    def sam_entries(self):
        for contig, reads in self.reads.items():
            ref_seq = self.genome[contig]
            for idx, sequence in enumerate(reads):
                yield '\t'.join([str(_) for _ in (
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
                    ''.join([self._phred(_) for _ in sequence.phred]),
                )])

    def save_to_sam(self, bam_filename, genome_filename):
        with NamedTemporaryFile(
                delete=False,
                mode='w',
                dir='.',
                suffix='.sam',
        ) as stream:
            sam_filename = stream.name
            stream.write(self.header())
            stream.write('\n')
            stream.write('\n'.join(self.sam_entries()))
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

    def _covered_seqs(self, contig_name, position):
        return [
            idx for idx, seq in enumerate(self[contig_name])
            if seq.start <= position < seq.stop
        ]

    @classmethod
    def _phred(cls, int_value):
        return chr(33 + int_value)

def ntf(*args, **kwargs):
    with NamedTemporaryFile(
            *args,
            delete=False,
            mode='w',
            **kwargs,
    ) as stream:
        filename = stream.name
    return filename
