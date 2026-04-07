"""Commandline tool for REDItools."""

import sys
from reditools.tools import analyze, annotate, find_repeats, index


def usage():
    """Print program usage."""
    print("""usage: reditools {analyze,find-repeats,index,annotate}

REDItools3

Run Modes:
  analyze            Find editing events in one or more alignment files.

  find-repeats       Find repetitive elements in a genome.

  index              Calculate editing indices from the output of `analyze`
                     mode.

  annotate           Annotate REDItools RNA output with DNA output
""")


if __name__ == '__main__':
    toolkit = {
        'analyze': analyze,
        'find-repeat': find_repeats,
        'index': index,
        'annotate': annotate,
    }
    if len(sys.argv) > 1:
        command = sys.argv.pop(1)
        if command in toolkit:
            toolkit[command].main()
        else:
            usage()
    else:
        usage()
