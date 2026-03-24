import argparse


def parse_options():  # noqa:WPS213,WPS432
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
        help='The bam file(s) to be analyzed.',
    )
    parser.add_argument(
        '-r',
        '--reference',
        help='Reference FASTA file.',
    )
    parser.add_argument(
        '-o',
        '--output-file',
        help='Path to write output to.',
        default='/dev/stdout',
    )
    parser.add_argument(
        '-s',
        '--strand',
        choices=(0, 1, 2),
        type=int,
        default=0,
        help='Strand: this can be 0 (unstranded),' +
        '1 (second strand oriented) or ' +
        '2 (first strand oriented).',
    )
    parser.add_argument(
        '-a',
        '--append-file',
        action='store_true',
        help='Appends results to file (and creates if not existing).',
    )
    parser.add_argument(
        '-g',
        '--region',
        help='Only analyzes the specified region.',
    )
    parser.add_argument(
        '-m',
        '--load-omopolymeric-file',
        help='BED file of omopolymeric positions.',
    )
    parser.add_argument(
        '-sf',
        '--splicing-file',
        help='The file containing splicing site positions.',
    )
    parser.add_argument(
        '-ss',
        '--splicing-span',
        type=int,
        default=4,
        help='The splicing span.',
    )
    parser.add_argument(
        '-mrl',
        '--min-read-length',
        type=int,
        default=30,
        help='Reads with length below -mrl will be discarded.',
    )
    parser.add_argument(
        '-q',
        '--min-read-quality',
        type=int,
        default=20,
        help='Reads with mapping quality below -q will be discarded.',
    )
    parser.add_argument(
        '-bq',
        '--min-base-quality',
        type=int,
        default=30,
        help='Base quality below -bq will bed discarded.',
    )
    parser.add_argument(
        '-mbp',
        '--min-base-position',
        type=int,
        default=0,
        help='Ignores the first -mbp bases in each read.',
    )
    parser.add_argument(
        '-Mbp',
        '--max-base-position',
        type=int,
        default=0,
        help='Ignores the last -Mpb bases in each read.',
    )
    parser.add_argument(
        '-l',
        '--min-read-depth',
        type=int,
        default=1,
        help='Only report on positions with at least -l read depth',
    )
    parser.add_argument(
        '-e',
        '--exclude-multis',
        default=False,
        help='Do not report any position with more than one alternate base.',
        action='store_true',
    )
    parser.add_argument(
        '-men',
        '--min-edits-per-nucleotide',
        type=int,
        default=0,
        help='Positions with fewer than -men edits will not be discarded.',
    )
    parser.add_argument(
        '-me',
        '--min-edits',
        type=int,
        default=1,
        help='The minimum number of editing events (per position). ' +
        'Positions with fewer than -me edits will be discarded.',
    )
    parser.add_argument(
        '-Men',
        '--max-editing-nucleotides',
        type=int,
        default=4,
        help='The maximum number of editing nucleotides, from 0 to 3 ' +
        '(per position). Positions whose columns have more than ' +
        '"max-editing-nucleotides" will not be included in the analysis.',
    )
    parser.add_argument(
        '-T',
        '--strand-confidence-threshold',
        type=float,
        default=0.7,
        help='Only report the strandedness if at least -T proportion of ' +
        'reads are of a given strand.',
    )
    parser.add_argument(
        '-C',
        '--strand-correction',
        default=False,
        help='Strand correction. Once the strand has been inferred, ' +
        'only bases according to this strand will be selected.',
        action='store_true',
    )
    parser.add_argument(
        '-V',
        '--verbose',
        default=False,
        help='Run in verbose mode.',
        action='store_true',
    )
    parser.add_argument(
        '-N',
        '--dna',
        default=False,
        help='Run REDItools on DNA-Seq data.',
        action='store_true',
    )
    parser.add_argument(
        '-B',
        '--bed_file',
        nargs='+',
        help='Only analyze regions in the provided BED file.',
    )
    parser.add_argument(
        '-t',
        '--threads',
        help='Number of threads for parallel processing.',
        type=int,
        default=1,
    )
    parser.add_argument(
        '--temp-dir',
        help='Location to save temporary files')
    parser.add_argument(
        '-w',
        '--window',
        help='How many bp should be processed by each thread at a time. ' +
        'Zero uses the full contig.',
        type=int,
        default=0,
    )
    parser.add_argument(
        '-k',
        '--exclude_regions',
        nargs='+',
        help='Skip regions in the provided BED file(s).',
    )
    parser.add_argument(
        '-E',
        '--exclude_reads',
        help='Text file listing read names to exclude from analysis.',
    )
    parser.add_argument(
        '-d',
        '--debug',
        default=False,
        help='Run in debug mode.',
        action='store_true',
    )
    parser.add_argument(
        '-v',
        '--variants',
        nargs='*',
        default=['all'],
        help='Which editing events to report. Edits should be two characters, '
        'separated by spaces. Use "all" to report all variants.',
    )

    return parser.parse_args()
