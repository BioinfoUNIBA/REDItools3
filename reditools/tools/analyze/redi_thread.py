import argparse

from reditools.region import Region
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.setup_alignment_manager import \
    setup_alignment_manager
from reditools.tools.analyze.setup_rtools import setup_rtools
from reditools.tools.analyze.write_results import write_results


class REDIThread:
    @classmethod
    def init(cls, options: argparse.Namespace) -> None:
        """Worker thread function for parallel REDItools analysis.

        Parameters
        ----------
        options : argparse.Namespace
            The command-line options.
        """
        cls.rtools = setup_rtools(options)
        cls.sam_manager = setup_alignment_manager(
            options.file,
            options.min_read_quality,
            options.min_read_length,
            options.exclude_reads,
        )
        cls.rtqc = RTChecks(options)
        cls.temp_dir = options.temp_dir

    @classmethod
    def analyze(
            cls,
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
        rtresults = cls.rtools.analyze(cls.sam_manager, region)
        return write_results(
            rtresults,
            cls.temp_dir,
            cls.rtqc,
            cls.rtools.log,
        )
