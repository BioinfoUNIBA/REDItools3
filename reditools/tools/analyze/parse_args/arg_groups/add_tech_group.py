"""CLI arguments for technical options."""
from __future__ import annotations

import tempfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse

def add_tech_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for technical options.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    tech_group = parser.add_argument_group(
        title="Technical Options",
    )
    tech_group.add_argument(
        "-V",
        "--verbose",
        default=False,
        help="Run in verbose mode.",
        action="store_true",
    )
    tech_group.add_argument(
        "-d",
        "--debug",
        default=False,
        help=(
            "Run in debug mode. Every step of REDItools logic will be printed "
            "to STDERR. If REDItools crashes, debug mode will print the stack "
            "trace as well."
        ),
        action="store_true",
    )
    tech_group.add_argument(
        "--temp-dir",
        help="Location to save temporary files",
        default=tempfile.gettempdir(),
    )
