"""Fast Logging for REDItools."""
import os
import socket
import sys
from datetime import datetime


class Logger(object):
    """Fast logger for REDItools."""

    silent_level = 'SILENT'
    info_level = 'INFO'
    debug_level = 'DEBUG'

    def __init__(self, level):
        """
        Create a new Logger.

        Parameters:
            level (str): either 'INFO' or 'DEBUG'
        """
        hostname = socket.gethostname()
        ip_addr = socket.gethostbyname(hostname)
        pid = os.getpid()
        self.hostname_string = f'{hostname}|{ip_addr}|{pid}'

        if level.upper() == self.debug_level:
            self.log = self._log_all
        elif level.upper() == self.info_level:
            self.log = self._log_info
        else:
            self.log = lambda *_: None

    def _log_all(self, level, message, *args):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        message = message.format(*args)
        sys.stderr.write(
            f'{timestamp} [{self.hostname_string}] ' +
            f'[{level}] {message}\n',
        )

    def _log_info(self, level, message, *args):
        if level == self.info_level:
            self._log_all(level, message, *args)
