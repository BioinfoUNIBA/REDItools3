"""CLI arguments for strand reporting and filtering."""
from __future__ import annotations

from typing import TYPE_CHECKING

from reditools import reditools
from reditools.tools.analyze.parse_args.bounded_types import bounded_float

if TYPE_CHECKING:
    import argparse

def add_strand_group(parser: argparse.ArgumentParser) -> None:
    """Add arguments for strand reporting and filtering.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI argparser
    """
    strand_group = parser.add_argument_group(
        title="Strandedness Options",
    )
    strand_group.add_argument(
        "-s",
        "--strand",
        choices=(
            reditools.UNSTRANDED_MODE,
            reditools.FORWARD_STRAND_MODE,
            reditools.REVERSE_STRAND_MODE,
        ),
        type=int,
        default=reditools.UNSTRANDED_MODE,
        help=(
            f"Infer RNA strand and filter reads not of the same strand. "
            f"This option may be {reditools.UNSTRANDED_MODE} (unstranded), "
            f"{reditools.FORWARD_STRAND_MODE} (read1 is original RNA), or "
            f"{reditools.REVERSE_STRAND_MODE} (read2 is original RNA). "
            "From RSeQC infer_experiment.py, 1++,1--,2+-,2-+ should be run "
            f"as --strand {reditools.FORWARD_STRAND_MODE} and 1+-,1-+,2++,2-- "
            f"should be run as --strand {reditools.REVERSE_STRAND_MODE}. "
            "From Salmon, forward libraries (ISF, MSF, OSF) should be run as "
            f"--strand {reditools.FORWARD_STRAND_MODE} and reverse libraries "
            f"(ISR, MSR, OSR) as --strand {reditools.REVERSE_STRAND_MODE}. "
            "All DNA sequencing experiments and non-stranded experiments "
            f"should be run with --strand {reditools.UNSTRANDED_MODE}."
        ),
    )
    strand_group.add_argument(
        "-T",
        "--strand-confidence-threshold",
        type=bounded_float(max_value=1),
        default=0.7,
        help=(
            "Only report the strandedness if at least -T proportion of "
            "reads are of a given strand. This option is only applicable "
            f"if -s/--strand is not {reditools.UNSTRANDED_MODE} and "
            "--separate-strands is not used."
        ),
    )
    strand_group.add_argument(
        "-C",
        "--strand-correction",
        default=False,
        help=(
            "Report the base complements for the Reference, AllSubs, and "
            "BaseCount columns in the output if the detected edit is on "
            "the minus strand. "
            "This option is only applicable if -s/--strand is not "
            f"{reditools.UNSTRANDED_MODE}."
        ),
        action="store_true",
    )
    strand_group.add_argument(
        "-p",
        "--separate-strands",
        default=False,
        help=(
            "Report edits for each strand on separate lines of the output. "
            "When using this option, all result filters will be applied to "
            "each strand individually. For example, a base position with 7 "
            "forward reads and 3 reverse reads would produce no output "
            " with --min-read-depth 10 and only report the fowrard reads "
            " with --min-read-depth 5. "
            "This option is only applicable if -s/--strand is not "
            f"{reditools.UNSTRANDED_MODE}."
        ),
        action="store_true",
    )
