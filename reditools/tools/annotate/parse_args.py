import argparse

__all__ = ('parse_args')


def parse_args():
    """
    Parse commandline options for REDItools.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(
        prog='reditools annotate',
        description='Annotates RNA REDItools output with DNA output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'rna_file',
        help='The REDItools output from RNA data',
    )
    parser.add_argument(
        'dna_file',
        help='The REDItools output from corresponding DNA data',
    )
    parser.add_argument(
        '-d',
        '--debug',
        help='Report stack trace on crash.',
        action='store_true',
    )
    order_group = parser.add_argument_group(
        title='Contig order options',
        description=(
            'By default, the annotate tool will determine the contig order '
            'by reading through the rna_file once before performing the '
            'annotation. Depending on the size of the rna_file, this could '
            'take a while. To skip this time consuming step, you can provide '
            'the contig order from another source.'
        ),
    )
    order_group.add_argument(
        '-b',
        '--bam',
        help='BAM file to get contig order from.',
    )
    order_group.add_argument(
        '-f',
        '--fai',
        help='FASTA Index file to get contig order from.',
    )
    options = parser.parse_args()

    if options.bam and options.fai:
        parser.error(
            message='Options -b/--bam and -f/--fai are mutually exclusive.',
        )

    return options
