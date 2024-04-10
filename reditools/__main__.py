"""Commandline tool for REDItools."""

import sys

from reditools import analyze, homopolymerics, index


def usage():
    """Print program usage."""
    print("""usage: reditools {analyze,find-repeats,index}

REDItools3

Run Modes:
  analyze            Find editing events in one or more alignment files.

  find-repeats       Find repetitive elements in a genome.

  index              Calculate editing indices from the output of `analyze`
                     mode.
""")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        command = sys.argv.pop(1)
        match command:
            case 'analyze':
                analyze.main()
            case 'find-repeats':
                homopolymerics.main()
            case 'index':
                index.main()
            case _:
                usage()
    else:
        usage()
