"""CLI argument parsing."""
from __future__ import annotations

import argparse
import json

from reditools import reditools
from reditools.tools.analyze.parse_args import arg_groups
from reditools.tools.analyze.parse_args.json_args import args_from_json


class DNAStrandError(argparse.ArgumentTypeError):
    """DNA mode requires strand set to 0."""

    def __init__(self) -> None:
        """Initialize self."""
        self.message = (
            "-N/--dna can only be used with "
            f"-s/--strand {reditools.UNSTRANDED_MODE}"
        )
        super().__init__(self.message)

class StrictConflictError(argparse.ArgumentTypeError):
    """Strict mode requires min-edits set to 1."""

    def __init__(self) -> None:
        """Initialize self."""
        self.message = "-S/--strict can only be used with -me/--min-edits 1."
        super().__init__(self.message)

def build_argument_parser() -> argparse.ArgumentParser:  # noqa: WPS213
    """Build the argument parser for reditools analyze.

    Returns
    -------
    argparse.ArgumentParser
        The configured argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="reditools analyze",
        description="REDItools3",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "file",
        nargs="+",
        help=(
            "The BAM file(s) to be analyzed. BAM files must be sorted and "
            "indexed."
        ),
    )
    parser.add_argument(
        "-r",
        "--reference",
        help=(
            "Reference genome FASTA file. (Note: REDItools runs fastest when "
            "BAM files have MD tags and -r is *not* used)."
        ),
    )
    parser.add_argument(
        "-g",
        "--region",
        help=(
            "Only analyzes the specified samtools formatted region. "
            "(1-index, start and end inclusive)."
        ),
    )
    parser.add_argument(
        "--resume",
        help=(
            "If REDItools crashed, will attempt to resume a stopped job from "
            "existing temporary files. Note: --temp-dir is required."
        ),
        action="store_true",
    )

    arg_groups.add_output_group(parser)
    arg_groups.add_bqf_group(parser)
    arg_groups.add_gr_group(parser)
    arg_groups.add_rf_group(parser)
    arg_groups.add_strand_group(parser)
    arg_groups.add_dna_group(parser)
    arg_groups.add_para_group(parser)
    arg_groups.add_tech_group(parser)
    arg_groups.add_leg_group(parser)

    return parser

def _delattr(namespace: argparse.Namespace, attribute: str) -> None:
    """Delete an attribute from a namespace.

    Parameters
    ----------
    namespace : argparse.Namespace
        Namespace to modify
    attribute : str
        Attribute to delete
    """
    delattr(namespace, attribute)  # noqa: WPS421

def fix_legacy_options(args: argparse.Namespace) -> None:
    """Adjust arguments to handle legacy options.

    The specific fixes are:
    - --dna sets --strand reditools.UNSTRANDED_MODE
    - --exclude-multis sets --max-editing-nucleotides 1
    - --strict sets --min-edits 1
    - --load-omopolymeric-file adds the input to --exclude-regions

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.

    Raises
    ------
    DNAStrandError, StrictConflictError
        If mutually exclusive options are provided.
    """
    if args.strand != reditools.UNSTRANDED_MODE and args.dna:
        raise DNAStrandError
    _delattr(args, "dna")

    if args.exclude_multis:
        args.max_editing_nucleotides = 1
    _delattr(args, "exclude_multis")

    if args.strict and args.min_edits != 1:
        raise StrictConflictError
    _delattr(args, "strict")

    if args.load_omopolymeric_file:
        if args.exclude_regions is None:
            args.exclude_regions = []
        args.exclude_regions.append(args.load_omopolymeric_file)
    _delattr(args, "load_omopolymeric_file")

def parse_args(sys_args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for reditools analyze.

    Parameters
    ----------
    sys_args : list[str] | None, optional
        The list of command-line arguments, by default None.

    Returns
    -------
    argparse.Namespace
        The parsed and validated command-line arguments.
    """
    parser = build_argument_parser()
    args = parser.parse_args(sys_args)

    if args.resume:
        temp_dir = args.temp_dir
        try:
            args = args_from_json(temp_dir)
        except (json.JSONDecodeError, OSError) as exc:
            parser.error(f"Unable to resume analysis.\n{exc}")
        args.resume = True
        args.temp_dir = temp_dir
        return args

    try:
        fix_legacy_options(args)
    except argparse.ArgumentTypeError as exc:
        parser.error(message=str(exc))

    if args.max_editing_nucleotides < args.min_edits:
        parser.error(
            "-Men/--max-editing-nucleotides cannot be smaller than "
            "-me/--min-edits.",
        )

    if args.strand == reditools.UNSTRANDED_MODE and args.strand_correction:
        parser.error(
            f"-s/--strand {reditools.UNSTRANDED_MODE} and "
            "-C/--strand-correction are mutually exclusive.",
        )

    if args.strand == 0 and args.separate_strands:
        parser.error(
            f"-s/--strand {reditools.UNSTRANDED_MODE} and "
            "-p/--separate-strands are mutually exclusive.",
        )

    return args

def args_to_string(args: argparse.Namespace) -> str:
    """Convert argparse options to a comma-separated string of key:value pairs.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command line options.

    Returns
    -------
    str
        A string representation of the options.
    """
    return ", ".join(
        [f"{_}:{getattr(args, _)}" for _ in vars(args)],  # noqa: WPS421
    )
