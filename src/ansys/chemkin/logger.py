# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Provides the singleton helper class for the logger."""

# logger from https://gist.github.com/huklee/cea20761dd05da7c39120084f52fcc7c
import datetime
import logging
from pathlib import Path
from typing import Dict


class SingletonType(type):
    """:no-index: Provides the singleton helper class for the logger."""

    _instances: Dict[type, type] = {}

    def __call__(cls, *args, **kwargs):
        """Call to redirect new instances to the singleton instance."""
        if cls not in cls._instances:
            cls._instances[cls] = super(SingletonType, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ChemkinLogger(object, metaclass=SingletonType):
    """Provides the singleton logger for the PyChemkin.

    Parameters
    ----------
    to_file : bool, default: False
        Whether to include the logs in a file.

    """

    def __init__(self, level: int = logging.ERROR, logger_name: str = "PyChemkin"):
        """Initialize the logger."""
        self._logger = logging.getLogger(logger_name)
        self._logger.setLevel(level)
        self._formatter = logging.Formatter(
            "%(asctime)s \t [%(levelname)s | %(filename)s:%(lineno)s] > %(message)s"
        )

    def get_logger(self) -> logging.Logger:
        """Get the logger.

        Returns
        -------
        Logger
            Logger.

        """
        return self._logger

    def set_level(self, level: int) -> None:
        """Set the logger output level.

        Parameters
        ----------
        level : int, {0, 10, 20, 30, 40, 50}
            Output Level of the logger.
            0 = NOTSET
            10 = DEBUG
            20 = INFO
            30 = WARNING
            40 = ERROR
            50 = CRITICAL

        """
        self._logger.setLevel(level=level)

    def enable_output(self, stream=None):
        """Enable logger output to a given stream.

        If a stream is not specified, ``sys.stderr`` is used.

        Parameters
        ----------
        stream: TextIO, default: ``sys.stderr``
            Stream to output the log output to.

        """
        # stdout
        stream_handler = logging.StreamHandler(stream)
        stream_handler.setFormatter(self._formatter)
        self._logger.addHandler(stream_handler)

    def add_file_handler(self, logs_dir: str = "./.log"):
        """Save logs to a file in addition to printing them to the standard output.

        Parameters
        ----------
        logs_dir : str, default: ``"./.log"``
            Directory of the logs.

        """
        now = datetime.datetime.now()
        if not Path.is_dir(Path(logs_dir)):
            Path.mkdir(Path(logs_dir))
        file_handler = logging.FileHandler(
            logs_dir + "/log_" + now.strftime("%Y-%m-%d") + ".log"
        )
        file_handler.setFormatter(self._formatter)
        self._logger.addHandler(file_handler)


logger = ChemkinLogger().get_logger()
ChemkinLogger().set_level(level=10)
ChemkinLogger().enable_output()
