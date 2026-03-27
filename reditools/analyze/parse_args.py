import argparse
import tempfile

__all__ = ('parse_args',)


def check_number_bounds(value, min=None, max=None):
    if min is not None and value < min:
        raise argparse.ArgumentTypeError(f'Value must be at least {min}.')
    if max is not None and value > max:
        raise argparse.ArgumentTypeError(
            f'Value cannot be larger than {max}.',
        )


def bounded_int(min=None, max=None):
    def subfn(value):
        try:
            value = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'invalid int value: {value}')
        check_number_bounds(value, min, max)
        return value
    return subfn


def bounded_float(min=None, max=None):
    def subfn(value):
        try:
            value = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'invalid float value: {value}')
        check_number_bounds(value, min, max)
        return value
    return subfn


def build_argument_parser():
    """
    Parse commandline options for REDItools.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(
        prog="reditools analyze",
        description='REDItools3',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'file',
        nargs='+',
        help=(
            'The BAM file(s) to be analyzed. BAM files must be sorted and '
            'indexed.'
        ),
    )
    parser.add_argument(
        '-r',
        '--reference',
        help=(
            'Reference genome FASTA file. (Note: REDItools runs fastest when '
            'BAM files have MD tags and -r is *not* used).'
        ),
    )
    parser.add_argument(
        '-g',
        '--region',
        help=(
            'Only analyzes the specified samtools formatted region. '
            '(1-index, start and end inclusive).'
        ),
    )
    output_group = parser.add_argument_group(
        title='Output Options',
    )
    output_group.add_argument(
        '-o',
        '--output-file',
        help='Path to write output to.',
        default='/dev/stdout',
    )
    output_group.add_argument(
        '-a',
        '--append-file',
        action='store_true',
        help='Appends results to file (and creates if not existing).',
    )
    bqf_group = parser.add_argument_group(
        title='Base/Read Quality Controls',
    )
    bqf_group.add_argument(
        '-mrl',
        '--min-read-length',
        type=int,
        default=30,
        help='Reads shorter than -mrl will be discarded.',
    )
    bqf_group.add_argument(
        '-q',
        '--min-read-quality',
        type=int,
        default=20,
        help='Reads with MAPQ below -q will be discarded.',
    )
    bqf_group.add_argument(
        '-bq',
        '--min-base-quality',
        type=int,
        default=30,
        help='Bases with a Phred quality score below -bq will bed discarded.',
    )
    bqf_group.add_argument(
        '-mbp',
        '--min-base-position',
        type=bounded_int(min=0),
        default=0,
        help='Ignores the first -mbp bases in each read.',
    )
    bqf_group.add_argument(
        '-Mbp',
        '--max-base-position',
        type=bounded_int(min=0),
        default=0,
        help='Ignores the last -Mpb bases in each read.',
    )
    bqf_group.add_argument(
        '-E',
        '--exclude-reads',
        help='Text file listing read names to exclude from analysis.',
    )
    bqf_group.add_argument(
        '--exclude_reads',
        help=argparse.SUPPRESS,
    )
    gr_group = parser.add_argument_group(
        title='Genomic Region Filters',
    )
    gr_group.add_argument(
        '-k',
        '--exclude-regions',
        nargs='+',
        help='Do not report on regions in the provided BED file(s).',
    )
    gr_group.add_argument(
        '--exclude_regions',
        nargs='+',
        help=argparse.SUPPRESS,
    )
    gr_group.add_argument(
        '-B',
        '--bed-file',
        nargs='+',
        help='Only reports on regions in the provided BED file(s).',
    )
    gr_group.add_argument(
        '--bed_file',
        nargs='+',
        help=argparse.SUPPRESS,
    )
    rf_group = parser.add_argument_group(
        title='Result Filters',
    )
    rf_group.add_argument(
        '-men',
        '--min-edits-per-nucleotide',
        type=int,
        default=0,
        help=(
            'Position where any variant has a frequency less than -men but '
            'more than zero will not be reported. (Corresponds to the '
            'BaseCount column.)'
        ),
    )
    rf_group.add_argument(
        '-me',
        '--min-edits',
        type=bounded_int(0, 4),
        default=1,
        help=(
            'Positions with fewer than -me unique variants (listed in the '
            'AllSubs column) will be excluded from the results.'
        ),
    )
    rf_group.add_argument(
        '-Men',
        '--max-editing-nucleotides',
        type=bounded_int(min=0, max=4),
        default=4,
        help=(
            'Positions with more than -Men unique variants (listed in the '
            'AllSubs column) will be excluded from the results.'
        ),
    )
    rf_group.add_argument(
        '-v',
        '--variants',
        nargs='*',
        default=['all'],
        help=(
            'Which editing events to report. Each edit should be two '
            'characters and separated by spaces (e.g. AG CT). Use "all" to '
            'report all variants. (Corresponds to the AllSubs column)'
        ),
    )
    rf_group.add_argument(
        '-l',
        '--min-read-depth',
        type=bounded_int(min=1),
        default=1,
        help=(
            'Only report on positions with at least -l reads (corresponds to '
            'the Coverage column.)'
        ),
    )

    strand_group = parser.add_argument_group(
        title='Strandedness Options',
    )
    strand_group.add_argument(
        '-s',
        '--strand',
        choices=(0, 1, 2),
        type=int,
        default=0,
        help=(
            'Strand can be 0 (unstranded), 1 (read1 is original RNA), or '
            '2 (read2 is original RNA). '
            'From RSeQC infer_experiment.py, 1++,1--,2+-,2-+ should be run '
            'as --strand 1 and 1+-,1-+,2++,2-- should be run as --strand 2. '
            'From Salmon, forward libraries (ISF, MSF, OSF) should be run as '
            '--strand 1 and reverse libraries (ISR, MSR, OSR) as --strand 2. '
            'All DNA sequencing experiments, single-end experiments, and '
            'nonstranded experiments should be run with --strand 0.'
        ),
    )
    strand_group.add_argument(
        '-T',
        '--strand-confidence-threshold',
        type=bounded_float(max=1),
        default=0.7,
        help=(
            'Only report the strandedness if at least -T proportion of '
            'reads are of a given strand. This option is only applicable '
            'if -s/--strand is not zero.'
        ),
    )
    strand_group.add_argument(
        '-C',
        '--strand-correction',
        default=False,
        help=(
            'Once the strand has been inferred, only bases according to this '
            'strand will be reported. This option is only applicable if '
            '-s/--strand is not zero.'
        ),
        action='store_true',
    )
    para_group = parser.add_argument_group(
        title='Parallel Processing Options',
    )
    para_group.add_argument(
        '-t',
        '--threads',
        help=(
            'Number of threads for parallel processing. Note that the '
            'maximum number of usable threads is equivalent to the number of '
            'chromosomes in your alignment genome unless you use the --window '
            'option.'
        ),
        type=bounded_int(min=1),
        default=1,
    )
    para_group.add_argument(
        '-w',
        '--window',
        help=(
            'How many bp should be processed by each thread at a time. '
            'Zero uses the full contig.'
        ),
        type=bounded_int(min=0),
        default=0,
    )
    tech_group = parser.add_argument_group(
        title='Technical Options',
    )
    tech_group.add_argument(
        '-V',
        '--verbose',
        default=False,
        help='Run in verbose mode.',
        action='store_true',
    )
    tech_group.add_argument(
        '-d',
        '--debug',
        default=False,
        help=(
            'Run in debug mode. Every step of REDItools logic will be printed '
            'to STDERR. If REDItools crashes, debug mode will print the stack '
            'trace as well.'
        ),
        action='store_true',
    )
    tech_group.add_argument(
        '--temp-dir',
        help='Location to save temporary files',
        default=tempfile.gettempdir(),
    )
    leg_group = parser.add_argument_group(
        title='Legacy Options',
    )
    leg_group.add_argument(
        '-e',
        '--exclude-multis',
        default=False,
        help=(
            'Do not report any position with more than one alternate base. '
            '(Equivalent to -Men/--max-editing-nucleotides 1)'
        ),
        action='store_true',
    )
    leg_group.add_argument(
        '-N',
        '--dna',
        default=False,
        help='Run REDItools on DNA-Seq data. (Equivalent to -s/--strand 0)',
        action='store_true',
    )
    leg_group.add_argument(
        '-m',
        '--load-omopolymeric-file',
        help=(
            'BED file of homopolymeric positions. Regions in the BED file '
            'will be excluded from the analysis. (Same effect as providing '
            'the BED file with the -k/--exclude-regions option.'
        ),
    )
    leg_group.add_argument(
        '-S',
        '--strict',
        help=(
            'Activate strict mode: only sites with edits will be included in '
            'the output. (Equivalent to -me/--min-edits 1)'
        ),
        action='store_true',
    )
    leg_group.add_argument(
        '-sf',
        '--splicing-file',
        help=(
            'The splicing file is a space delimited file with five columns: '
            'chromosome, start position (one-index inclusive), stop '
            '(ignored),  splice (either A for acceptor or D for donor), and '
            'strand (either + or -). A header is optional, but must start '
            'with #. Used in conjunctions with -ss/--splicing-span.'
        ),
    )
    leg_group.add_argument(
        '-ss',
        '--splicing-span',
        type=bounded_int(min=1),
        default=4,
        help=(
            'The splicing span. Used in conjunction with -sf/--splicing-file.'
        ),
    )

    return parser


def check_dna_mode(args):
    if args.strand != 0 and args.dna:
        raise Exception('-N/--dna can only be used with -s/--strand 0.')
    delattr(args, 'dna')


def check_exclude_multis(args):
    if args.exclude_multis:
        setattr(args, 'max_editing_nucleotides', 1)
    delattr(args, 'exclude_multis')


def check_strict_mode(args):
    if args.strict:
        if args.min_edits != 1:
            raise Exception(
                '-S/--strict can only be used with -me/--min-edits 1.'
            )
    delattr(args, 'strict')


def check_load_omopolymeric_file(args):
    if args.load_omopolymeric_file:
        if args.exclude_regions is None:
            setattr(args, 'exclude_regions', [args.load_omopolymeric_file])
        else:
            args.exclude_regions.append(args.load_omopolymeric_file)
    delattr(args, 'load_omopolymeric_file')


def test_edit_frequency(args):
    if args.max_editing_nucleotides < args.min_edits:
        raise Exception(
            '-Men/--max-editing-nucleotides cannot be smaller than '
            '-me/--min-edits.',
        )


def test_strand_args(args):
    if args.strand == 0 and args.strand_correction:
        raise Exception(
            '-s/--strand 0 and -C/--strand-correction are mutually exclusive.'
        )


def parse_args(sys_args=None):
    parser = build_argument_parser()
    args = parser.parse_args(sys_args)
    try:
        check_dna_mode(args)
        check_exclude_multis(args)
        check_strict_mode(args)
        check_load_omopolymeric_file(args)

        test_edit_frequency(args)
        test_strand_args(args)
    except Exception as e:
        parser.error(message=str(e))

    return args
