"""CLI entry point for REDItools."""

import sys

from reditools.tools.analyze import main as analyze
from reditools.tools.annotate import main as annotate
from reditools.tools.find_repeats import main as find_repeats
from reditools.tools.index import main as index


def usage() -> None:
    """Print the usage information for the REDItools3 toolkit."""
    usage_str = """usage: reditools {analyze,find-repeats,index,annotate}

REDItools3

Run Modes:
  analyze            Find editing events in one or more alignment files.

  find-repeats       Find repetitive elements in a genome.

  index              Calculate editing indices from the output of `analyze`
                     mode.

  annotate           Annotate REDItools RNA output with DNA output
"""
    print(usage_str)  # noqa: WPS421 T201


if __name__ == "__main__":
    toolkit = {
        "analyze": analyze,
        "find-repeats": find_repeats,
        "index": index,
        "annotate": annotate,
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
