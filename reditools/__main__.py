"""Commandline tool for REDItools."""

import sys

from reditools.tools import analyze, annotate, find_repeats, index


def usage():
    """Print program usage."""
    usage_str = """usage: reditools {analyze,find-repeats,index,annotate}

REDItools3

Run Modes:
  analyze            Find editing events in one or more alignment files.

  find-repeats       Find repetitive elements in a genome.

  index              Calculate editing indices from the output of `analyze`
                     mode.

  annotate           Annotate REDItools RNA output with DNA output
"""
    print(usage_str)  # noqa: WPS421


if __name__ == '__main__':
    toolkit = {
        'analyze': analyze,
        'find-repeats': find_repeats,
        'index': index,
        'annotate': annotate,
    }
    if len(sys.argv) > 1:
        command = sys.argv.pop(1)
        tool = toolkit.get(command)
        if tool is None:
            usage()
        else:
            tool.main()
    else:
        usage()
