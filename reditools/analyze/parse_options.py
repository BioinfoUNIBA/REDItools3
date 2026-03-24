import argparse

__all__ = ('parse_options',)


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


def test_dna_strand_conflict(options):
    if options.strand != 0 and options.dna:
        raise Exception('Options --dna and --strand are mutually exclusive.')


def test_multis_conflict(options):
    if options.exclude_multis and options.max_editing_nucleotides != 1:
        raise Exception(
            'Options --exclude-multis and --max-editing-nucleotides are '
            'mutually exclusive.',
        )


def test_edit_frequency(options):
    if options.max_editing_nucleotides < options.min_edits:
        raise Exception(
            '-Men/--max-editing-nucleotides cannot be smaller than '
            '-me/--min-edits.',
        )


def parse_options():  # noqa:WPS213
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
        help='The BAM file(s) to be analyzed.',
    )
    parser.add_argument(
        '-r',
        '--reference',
        help=(
            'Reference genome FASTA file. (Note: REDItools runs fastest when '
            'BAM files have MD tags and -r is *not* used).'
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
        default=30,  # noqa:WPS432
        help='Reads shorter than -mrl will be discarded.',
    )
    bqf_group.add_argument(
        '-q',
        '--min-read-quality',
        type=int,
        default=20,  # noqa:WPS432
        help='Reads with MAPQ below -q will be discarded.',
    )
    bqf_group.add_argument(
        '-bq',
        '--min-base-quality',
        type=int,
        default=30,  # noqa:WPS432
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
        help='Skip regions in the provided BED file(s).',
    )
    gr_group.add_argument(
        '--exclude_regions',
        help=argparse.SUPPRESS,
    )
    gr_group.add_argument(
        '-g',
        '--region',
        help=(
            'Only analyzes the specified SAM region '
            '(1-index, start and end inclusive).'
        ),
    )
    gr_group.add_argument(
        '-B',
        '--bed-file',
        nargs='+',
        help='Only reports on regions in the provided BED file.',
    )
    gr_group.add_argument(
        '--bed_file',
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
        help='
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
        default=4,  # noqa:WPS432
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
            'report all variants.'
        ),
    )
    rf_group.add_argument(
        '-l',
        '--min-read-depth',
        type=bounded_int(min=1),
        default=1,
        help=(
            'Only report on positions with at least -l reads (corresponds to '
            'the Coverage column.) This is calculated after all other filters '
            'have been applied.'
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
            'Strand: this can be 0 (unstranded),'
            '1 (second strand oriented) or '
            '2 (first strand oriented).'
        ),
    )
    strand_group.add_argument(
        '-T',
        '--strand-confidence-threshold',
        type=bounded_float(max=1),
        default=0.7,  # noqa:WPS432
        help=(
            'Only report the strandedness if at least -T proportion of '
            'reads are of a given strand.'
        ),
    )
    strand_group.add_argument(
        '-C',
        '--strand-correction',
        default=False,
        help=(
            'Strand correction. Once the strand has been inferred, '
            'only bases according to this strand will be selected.'
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
        '-sf',
        '--splicing-file',
        help=(
            'The splicing file is a space delimited file with columns '
            'chromosome, start position (zero-index inclusive), splice '
            '(either D or A), and strand (either + or -). A header is '
            'optional, but must start with #.'
        ),
    )
    leg_group.add_argument(
        '-ss',
        '--splicing-span',
        type=bounded_int(min=1),
        default=4,
        help='The splicing span (used in conjunction with --splicing-file.)',
    ) 

    options = parser.parse_args()

    try:
        test_dna_strand_conflict(options)
        test_multis_conflict(options)
        test_edit_frequency(options)
    except Exception as e:
        parser.error(message=str(e))

    return options
