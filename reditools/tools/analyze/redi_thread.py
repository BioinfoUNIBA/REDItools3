import argparse

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
    ) -> str:
        """Analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.

        Returns
        -------
        str
            The path to the temporary file containing the results.
        """
        rtresults = self.rtools.analyze(self.sam_manager, region)
        return write_results(
            rtresults,
            self.temp_dir,
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
    def analyze(cls, region: Region) -> str:
        """Instruct thread to analyze a specific genomic region.

        Parameters
        ----------
        region : Region
            The genomic region to analyze.

        Returns
        -------
        str
            The path to the temporary file containing the results.
        """

        if cls.thread is None:
            raise AttributeError('REDIThreadManager not initialized.')
        return cls.thread.analyze(region)
