import argparse
from pathlib import Path

from reditools.region import Region
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.setup_alignment_manager import \
    setup_alignment_manager
from reditools.tools.analyze.setup_rtools import setup_rtools
from reditools.tools.analyze.write_results import write_results


class REDIThread:
    def __init__(self, options: argparse.Namespace) -> None:
        """Worker thread function for parallel REDItools analysis.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        self.rtools = setup_rtools(options)
        self.sam_manager = setup_alignment_manager(
            options.file,
            options.min_read_quality,
            options.min_read_length,
            options.exclude_reads,
        )
        self.rtqc = RTChecks(options)
        self.temp_dir = options.temp_dir

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
        return write_results(
            rtresults,
            filename,
            self.rtqc,
            self.rtools.log,
        )

class REDIThreadManager:
    """Manages a worker thread function for parallel REDItools analysis."""

    thread: None | REDIThread = None

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
    def analyze(cls, region: Region, filename: str) -> None:
        """Instruct thread to analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.
        filename : str
            Path to save output to.
        """

        if cls.thread is None:
            raise AttributeError('REDIThreadManager not initialized.')
        done_file = f'{filename}.done'
        if not Path(done_file).exists():
            cls.thread.analyze(region, filename)
            Path(done_file).touch()
