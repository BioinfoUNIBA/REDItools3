"""Fast Logging for REDItools."""
import os
import socket
import sys
from datetime import datetime
from typing import Any


class Logger:
    """Fast logger for REDItools."""

    silent_level = 'SILENT'
    info_level = 'INFO'
    debug_level = 'DEBUG'

    def __init__(self, level: str):
        """
        Create a new Logger.

        Parameters:
            level (str): either 'INFO' or 'DEBUG'
        """
        hostname = socket.gethostname()
        ip_addr = socket.gethostbyname(hostname)
        pid = os.getpid()
        self.hostname_string = f'{hostname}|{ip_addr}|{pid}'
        self._level = level.upper()

        if self._level == self.debug_level:
            self.log = self._log_all
        elif self._level == self.info_level:
            self.log = self._log_info
        else:
            self.log = self._log_silent

    @property
    def level(self):
        return self._level

    def _log_all(self, level: str, message: str, *args: Any):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        message = message.format(*args)
        sys.stderr.write(
            f'{timestamp} [{self.hostname_string}] ' +
            f'[{level}] {message}\n',
        )

    def _log_info(self, level: str, message: str, *args: Any):
        if level == self.info_level:
            self._log_all(level, message, *args)

    def _log_silent(self, level: str, message: str, *args: Any):
        pass  # noqa: WPS420
