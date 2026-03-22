from reditools.logger import logger

class QCCheck:
    @staticmethod
    def run_check(rtools, bases):
        raise NotImplementedError()

    @staticmethod
    def _log(rtools, msg, *args):
        rtools.log(
            Logger.debug_level,
            msg,
            *args,
        )
