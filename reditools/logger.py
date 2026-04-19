import os
import socket
import sys
from datetime import datetime
from typing import Any


class Logger:
    """
    Handle logging operations with different severity levels.

    Attriutes
    ----------
    silent_level : str = 'SILENT'
        Do not output anything
    info_level : str = 'INFO'
        Only output messages if the log level is info_level
    debug_level : str = 'DEBUG'
        Output all messages
    """

    silent_level = 'SILENT'
    info_level = 'INFO'
    debug_level = 'DEBUG'

    def __init__(self, level: str):
        """
        Initialize the Logger with a specified logging level.

        Parameters
        ----------
        level : str
            The logging level ('SILENT', 'INFO', or 'DEBUG').
        """
        hostname = socket.gethostname()
        ip_addr = socket.gethostbyname(hostname)
        pid = os.getpid()
        self.hostname_string = f'{hostname}|{ip_addr}|{pid}'
        self._level = level.upper()

        if self._level == self.debug_level:
            self._log_fn = self._log_all
        elif self._level == self.info_level:
            self._log_fn = self._log_info
        else:
            self._log_fn = self._log_silent

    def log(self, level: str, message: str, *args: Any):
        """Conditionally output a message to STDERR.

        Parameters
        ----------
        level : str
            The level for the message.
        message : str
            The message for output.
        *args : str
            Elements to fill in the message using fstring formatting.
        """
        self._log_fn(level, message, *args)

    @property
    def level(self) -> str:
        """
        Get the current logging level.

        Returns
        -------
        str
            The upper-case string representing the logging level.
        """
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
