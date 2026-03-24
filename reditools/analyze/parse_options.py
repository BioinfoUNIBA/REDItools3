import argparse
import sys


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
            'Reference FASTA file. (Note: REDItools runs fastest when BAM '
            'files have MD tags and -r is *not* used).'
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
        help='Bases Phred quality score below -bq will bed discarded.',
    )
    bqf_group.add_argument(
        '-mbp',
        '--min-base-position',
        type=int,
        default=0,
        help='Ignores the first -mbp bases in each read.',
    )
    bqf_group.add_argument(
        '-Mbp',
        '--max-base-position',
        type=int,
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
        '-m',
        '--load-omopolymeric-file',
        help='BED file of omopolymeric positions.',
    )
    gr_group.add_argument(
        '-g',
        '--region',
        help='Only analyzes the specified region.',
    )
    gr_group.add_argument(
        '-B',
        '--bed-file',
        nargs='+',
        help='Only analyze regions in the provided BED file.',
    )
    gr_group.add_argument(
        '--bed_file',
        help=argparse.SUPPRESS,
    )
    rf_group = parser.add_argument_group(
        title='Result Filters',
    )
    rf_group.add_argument(
        '-l',
        '--min-read-depth',
        type=int,
        default=1,
        help='Only report on positions with at least -l read depth',
    )
    rf_group.add_argument(
        '-men',
        '--min-edits-per-nucleotide',
        type=int,
        default=0,
        help='Positions with fewer than -men edits will be discarded.',
    )
    rf_group.add_argument(
        '-me',
        '--min-edits',
        type=int,
        default=1,
        help='The minimum number of editing events (per position). ' +
        'Positions with fewer than -me edits will be discarded.',
    )
    rf_group.add_argument(
        '-Men',
        '--max-editing-nucleotides',
        type=int,
        default=4,  # noqa:WPS432
        help='The maximum number of editing nucleotides, from 0 to 4 ' +
        '(per position). Positions whose columns have more than ' +
        '"max-editing-nucleotides" will not be included in the analysis.',
    )
    rf_group.add_argument(
        '-v',
        '--variants',
        nargs='*',
        default=['all'],
        help='Which editing events to report. Edits should be two characters, '
        'separated by spaces. Use "all" to report all variants.',
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
        help='Strand: this can be 0 (unstranded),' +
        '1 (second strand oriented) or ' +
        '2 (first strand oriented).',
    )
    strand_group.add_argument(
        '-T',
        '--strand-confidence-threshold',
        type=float,
        default=0.7,  # noqa:WPS432
        help='Only report the strandedness if at least -T proportion of ' +
        'reads are of a given strand.',
    )
    strand_group.add_argument(
        '-C',
        '--strand-correction',
        default=False,
        help='Strand correction. Once the strand has been inferred, ' +
        'only bases according to this strand will be selected.',
        action='store_true',
    )
    para_group = parser.add_argument_group(
        title='Parallel Processing Options',
    )
    para_group.add_argument(
        '-t',
        '--threads',
        help='Number of threads for parallel processing.',
        type=int,
        default=1,
    )
    para_group.add_argument(
        '-w',
        '--window',
        help='How many bp should be processed by each thread at a time. ' +
        'Zero uses the full contig.',
        type=int,
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
        help='Run in debug mode.',
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
            '(Equivalent to --max-editing-nucleotides 1)'
        ),
        action='store_true',
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
        type=int,
        default=4,
        help='The splicing span (used in conjunction with --splicing-file.)',
    )
    parser.add_argument(
        '-N',
        '--dna',
        default=False,
        help='Run REDItools on DNA-Seq data. (Equivalent to --strand 0)',
        action='store_true',
    )

    print(parser.parse_args())
    sys.exit(1)

    return parser.parse_args()
