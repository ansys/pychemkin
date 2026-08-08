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

"""Chemkin profile keyword data model."""

import copy

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import P_ATM
from ansys.chemkin.core.keyword import Keyword
from ansys.chemkin.core.logger import logger


class Profile:
    """Chemkin profile keyword class."""

    def __init__(
        self,
        key: str,
        x: npt.NDArray[np.double],
        y: npt.NDArray[np.double],
        label: bool = False,
    ):
        """Create a profile object."""
        # initialization
        self._profilekeyword = ""
        self._status = 0
        # check
        if key in Keyword.profilekeywords:
            self._profilekeyword = key
        elif label:
            this_key = key.split()[0]
            if this_key in Keyword.profilekeywords:
                self._profilekeyword = key
        else:
            msg = [
                Color.PURPLE,
                "profile is not available under the reactor model",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            self._status = -1
            return
        # profile data sizes
        xsize = len(x)
        ysize = len(y)
        if xsize == ysize:
            self._size = xsize
            # independent variable (time, location, grid, ...)
            if isinstance(x, np.double):
                self._pos = copy.deepcopy(x)
            else:
                self._pos = np.array(x, dtype=np.double)
            # dependent variable value at the corresponding position
            if isinstance(y, np.double):
                self._val = copy.deepcopy(y)
            else:
                self._val = np.array(y, dtype=np.double)
        else:
            msg = [
                Color.PURPLE,
                "the number of positions does not match the number of values\n",
                Color.SPACEx6,
                "number of positions =",
                str(xsize),
                "\n",
                Color.SPACEx6,
                "number of values    =",
                str(ysize),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            self._status = -2

    @property
    def size(self) -> int:
        """Get number of data points in the profile."""
        return self._size

    @property
    def status(self) -> int:
        """Get the validity of the profile object."""
        return self._status

    @property
    def pos(self) -> npt.NDArray[np.double]:
        """Get position values of profile data."""
        return self._pos

    @property
    def value(self) -> npt.NDArray[np.double]:
        """Get variable values of profile data."""
        return self._val

    @property
    def profilekey(self) -> str:
        """Get profile keyword."""
        return self._profilekeyword

    def show(self):
        """Show the profile data."""
        print(f"profile size: {self._size:d}")
        print(f" position           {self._profilekeyword:s}  ")
        for i in range(self._size):
            print(f"{self._pos[i]:f}         {self._val[i]}")

    def resetprofile(
        self, size: int, x: npt.NDArray[np.double], y: npt.NDArray[np.double]
    ):
        """Reset the profile data."""
        # check array size
        if size == self._size:
            # new profile has the same size
            self._pos[:] = x[:]
            self._val[:] = y[:]
        else:
            # new profile has different size
            self._size = size
            # resize the arrays
            self._pos.resize(size, refcheck=False)
            self._val.resize(size, refcheck=False)
        # fill the arrays with new values
        self._pos[:] = x[:]
        self._val[:] = y[:]

    def getprofile_as_string_list(self) -> tuple[int, list[str]]:
        """Create the keyword input lines as a list for Chemkin applications."""
        # initialization
        line = []
        # special treatment for pressure profile
        factor = 1.0e0
        if self._profilekeyword == "PPRO":
            if not Keyword.no_fullkeyword:
                # use Full Keywords: pressure units = atm
                factor = P_ATM
        # assemble profile keyword lines
        for i in range(self._size):
            thisline = (
                self._profilekeyword
                + Keyword.fourspaces
                + str(self._pos[i])
                + Keyword.fourspaces
                + str(self._val[i] / factor)
            )
            line.append(thisline)
        return self._size, line
