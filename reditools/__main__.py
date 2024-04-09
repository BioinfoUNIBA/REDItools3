"""Commandline tool for REDItools."""

import argparse
import sys

from reditools import analyze, homopolymerics, index


def usage():
    """Print program usage."""
    parser = argparse.ArgumentParser(
        description='REDItools 2.0',
        prog='reditools',
    )
    parser.add_argument(
        'command',
        choices=[
            'analyze',
            'find-repeats',
        ],
        help='Either analyze or find-repeats.',
    )
    parser.parse_args()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        command = sys.argv.pop(1)
        if command == 'analyze':
            analyze.main()
        elif command == 'find-repeats':
            homopolymerics.main()
        elif command == 'index':
            index.main()
        else:
            usage()
    else:
        usage()
