"""Create and manage threads for REDItools analyze tool."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pysam

from reditools import reditools
from reditools.logger import Logger
from reditools.rtannotater import RTAnnotater
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.setup_alignment_manager import (
    setup_alignment_manager,
)
from reditools.tools.analyze.write_results import write_results

if TYPE_CHECKING:
    import argparse

    from reditools.region import Region

class UninitializedError(AttributeError):
    """REDIThread.init has not been called."""

    def __init__(self) -> None:
        """Initialize self."""
        self.message = "REDIThreadManager not initialized."
        super().__init__(self.message)

def log_level(options: argparse.Namespace) -> str:
    """Determine Logger log level from CLI options.

    Parameters
    ----------
    options : argparse.Namespace
        CLI options

    Returns
    -------
    int : Logger log level
    """
    if options.debug:
        return Logger.debug_level
    if options.verbose:
        return Logger.info_level
    return Logger.silent_level

def contig_order_from_bam(bam_fname: str) -> dict[str, int]:
    """Get contig order from a BAM file.

    Parameters
    ----------
    bam_fname : str
        The path to the BAM file.

    Returns
    -------
    dict[str, int]
        A dictionary mapping contig names to their 1-based order.
    """
    with pysam.AlignmentFile(bam_fname, ignore_truncation=True) as bam:
        return {
            contig: idx
            for idx, contig in enumerate(bam.references, start=1)
        }


class REDIThread:
    """Worker thread for parallel REDItools analysis."""

    def __init__(self, options: argparse.Namespace) -> None:
        """Worker thread function for parallel REDItools analysis.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        self.rtools = reditools.REDItools()
        self.rtools.set_log_level(log_level(options))

        if options.reference:
            self.rtools.add_reference(options.reference)

        self.rtools.min_base_position = options.min_base_position
        self.rtools.max_base_position = options.max_base_position
        self.rtools.min_base_quality = options.min_base_quality

        self.rtools.strand = options.strand
        self.rtools.strand_confidence_threshold = \
            options.strand_confidence_threshold

        if options.strand_correction:
            self.rtools.use_strand_correction()

        self.sam_manager = setup_alignment_manager(
            options.file,
            options.min_read_quality,
            options.min_read_length,
            options.exclude_reads,
        )
        self.rtqc = RTChecks(options)
        self.temp_dir = options.temp_dir
        self.number_strand_output = options.number_strand_output

    def analyze(
            self,
            region: Region,
            filename: str,
    ) -> None:
        """Analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.
        filename : str
            Path to save output to.
        """
        rtresults = self.rtools.analyze(self.sam_manager, region)
        write_results(
            rtresults,
            filename,
            self.rtqc,
            self.rtools.log,
            self.number_strand_output,
        )

class REDIThreadDNA(REDIThread):
    """Worker thread for parallel REDItools analysis with both RNA and DNA."""

    def __init__(self, options: argparse.Namespace) -> None:
        """Worker thread function for parallel REDItools analysis.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        super().__init__(options)
        self.dna_rtools = reditools.REDItools()
        self.dna_rtools.set_log_level(log_level(options))

        if options.reference:
            self.dna_rtools.add_reference(options.reference)

        self.dna_rtools.min_base_position = options.dna_min_base_position
        self.dna_rtools.max_base_position = options.dna_max_base_position
        self.dna_rtools.min_base_quality = options.dna_min_base_quality
        self.dna_rtools.strand = reditools.UNSTRANDED_MODE

        self.dna_sam_manager = setup_alignment_manager(
            [options.dna_sample],
            options.dna_min_read_quality,
            options.dna_min_read_length,
            options.exclude_reads,
        )
        self.dna_rtqc = RTChecks(options)
        self.filter = options.dna_filter

        with pysam.AlignmentFile(
                options.dna_sample,
                ignore_truncation=True,
        ) as bam:
            contig_order = {
                contig: idx
                for idx, contig in enumerate(bam.references, start=1)
            }

        self.rta = RTAnnotater(contig_order)

    def analyze(
            self,
            region: Region,
            filename: str,
    ) -> None:
        """Analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.
        filename : str
            Path to save output to.
        """
        rna_results = self.rtools.analyze(self.sam_manager, region)
        dna_results = self.dna_rtools.analyze(self.dna_sam_manager, region)
        if self.filter:
            rtresults = self.rta.filtered_annotate(rna_results, dna_results)
        else:
            rtresults = self.rta.annotate(rna_results, dna_results)
        write_results(
            rtresults,
            filename,
            self.rtqc,
            self.rtools.log,
            self.number_strand_output,
        )


class REDIThreadManager:
    """Manages a worker thread function for parallel REDItools analysis."""

    thread: REDIThread | None = None

    @classmethod
    def init_thread(cls, options: argparse.Namespace) -> None:
        """Initialize a REDIThread.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        cls.thread = REDIThread(options)

    @classmethod
    def init_dna_thread(cls, options: argparse.Namespace) -> None:
        """Initialize a REDIThreadDNA.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        cls.thread = REDIThreadDNA(options)

    @classmethod
    def analyze(cls, region: Region, filename: str) -> bool:
        """Instruct thread to analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.
        filename : str
            Path to save output to.

        Returns
        -------
        bool
            True if the analysis runs, False if previously completed.
        """
        if cls.thread is None:
            raise UninitializedError
        done_file = f"{filename}.done"
        if Path(done_file).exists():
            return False
        cls.thread.analyze(region, filename)
        Path(done_file).touch()
        return True
