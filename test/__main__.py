import unittest
from .region import TestRegion
from .fasta_file import TestRTFastaFile
from .region_collection import TestRegionCollection
from .file_utils import TestFileUtils
from .compiled_position import TestCompiledPosition
from .analyze.region_args import TestRegionArgs
from .analyze.parse_args_utils import TestParseArgsUtils
from .analyze.parse_args import TestParseArgs
from .alignment_file import TestRTAlignmentFile
from .alignment_manager import TestAlignmentManager
from .compiled_reads import TestCompiledReads
from .rtchecks import TestRTChecks

unittest.main()
