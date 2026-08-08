# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

"""Chemkin keyword models and helpers."""

from typing import Union

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger


class Keyword:
    """A Chemkin style keyword."""

    # supported Chemkin keyword data types
    _keyworddatatypes = ["bool", "int", "float", "str"]
    _valuetypes = (bool, int, float, str)
    # required keywords that are given as reactor properties or as mixture properties
    # and will be set by using the KINAll0D_SetupBatchInputs call
    _protectedkeywords = [
        "CONP",
        "CONV",
        "TRAN",
        "STST",
        "TGIV",
        "ENRG",
        "PRES",
        "TEMP",
        "TAU",
        "TIME",
        "XEND",
        "FLRT",
        "VDOT",
        "SCCM",
        "DIAM",
        "REAC",
        "GAS",
        "INIT",
        "XEST",
        "SURF",
        "ACT",
        "TINL",
        "FUEL",
        "OXID",
        "PROD",
        "ASEN",
        "ATLS",
        "RTLS",
        "EPST",
        "EPSS",
    ]
    gasspecieskeywords = ["REAC", "XEST", "FUEL", "OXID"]
    flowratekeywords = ["FLRT", "VDOT", "VEL", "SCCM"]
    profilekeywords = [
        "TPRO",
        "PPRO",
        "VPRO",
        "QPRO",
        "AINT",
        "AEXT",
        "DPRO",
        "FPRO",
        "SCCMPRO",
        "VDOTPRO",
        "VELPRO",
        "TINPRO",
        "AFLO",
    ]
    fourspaces = "    "
    # Under the default API-call mode, important keywords (the _protectedkeywords)
    # are set by direct API calls, the rest of the keywords can be set by
    # keyword input lines (i.e., using the setkeyword method).
    # Under the full-keyword mode, all keywords and their parameters are set
    # by keyword input lines, and specifying those _protectedkeywords via
    # the setkeyword method are required.
    no_fullkeyword = True  # default: API-call mode

    def __init__(self, phrase: str, value: Union[float, bool, str], data_type: str):
        """Initialize the Chemkin keyword."""
        self._set = False
        ierr = 0
        # check value data type
        if data_type not in Keyword._keyworddatatypes:
            # the declared data type is not supported
            msg = [
                Color.PURPLE,
                "unsupported data type specified",
                data_type,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            if not isinstance(value, (bool, int, float, str)):
                # value does not match the declared data type
                msg = [
                    Color.PURPLE,
                    "variable has different data type",
                    str(type(value)),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            ierr = 1
        # block the protected keywords
        if Keyword.no_fullkeyword:
            if phrase.upper() in Keyword._protectedkeywords:
                msg = [
                    Color.PURPLE,
                    "use reactor property setter to assign",
                    phrase,
                    "value\n",
                    Color.SPACEx6,
                    "for example, to set the reactor volume use:",
                    '"MyBatchReactor.volume = 100"',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 2
        if ierr > 0:
            return
        self._key = phrase  # Chemkin keyword phrase
        self._value = value  # value assigned to the keyword
        self._data_type = data_type  # data type of the values
        self._prefix = ""  # a prefix to the keyword that can be used
        # to comment out/disable the keyword by setting it to '!'
        self._set = True

    @staticmethod
    def setfullkeywords(mode: bool):
        """Require all keywords and their parameters must be specified."""
        if mode:
            # turn ON the full keyword mode (no checking on protected keywords)
            Keyword.no_fullkeyword = False
        else:
            Keyword.no_fullkeyword = True

    def show(self):
        """Display the Chemkin keyword and its parameter value."""
        if self._set:
            if isinstance(self._value, (int, float)):
                msg = [
                    Color.YELLOW,
                    "keyword",
                    "'" + self._key + "':",
                    "value =",
                    str(self._value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            elif isinstance(self._value, bool):
                if self._value:
                    msg = [
                        Color.YELLOW,
                        "keyword",
                        "'" + self._key + "':",
                        "value = True",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                else:
                    msg = [
                        Color.YELLOW,
                        "keyword",
                        "'" + self._prefix + self._key + "':",
                        "value = Disabled",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
            else:
                msg = [
                    Color.YELLOW,
                    "keyword",
                    "'" + self._key + "':",
                    "value =",
                    str(self._value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        else:
            msg = [
                Color.YELLOW,
                "keyword",
                "'" + self._key + "':",
                "value not set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    def resetvalue(self, value: Union[float, bool, str]):
        """Reset the parameter value of an existing keyword."""
        if isinstance(value, Keyword._valuetypes):
            if isinstance(value, bool):
                if value:
                    # true: keep the keyword active
                    self._prefix = ""
                    self._value = True
                else:
                    # false: disable the keyword
                    self._prefix = "!"
                    self._value = False
            else:
                # integer, float, or string parameter
                self._value = value
        else:
            msg = [
                Color.PURPLE,
                "value has a wrong data type",
                type(value),
                "value will not be reset.\n",
                Color.SPACEx6,
                "data type expected by keyword",
                self._key,
                "is",
                self._data_type,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def parametertype(self) -> type:
        """Get parameter type of the keyword."""
        return type(self._value)

    @property
    def value(self) -> Union[int, float, bool, str]:
        """Get parameter value of the keyword."""
        # extract the keyword value
        if self._data_type == "bool":
            mode = self._prefix != "!"
            return mode
        return self._value

    @property
    def keyphrase(self) -> str:
        """Get the phrase of the keyword."""
        return self._key

    @property
    def keyprefix(self) -> bool:
        """Get the status of the keyword."""
        return self._prefix != "!"

    def getvalue_as_string(self) -> tuple[int, str]:
        """Create the keyword input line for Chemkin applications."""
        # initialization
        line = ""
        linelength = 0
        # assembly the keyword line
        if self._data_type == "bool":
            # boolean keyword (active or disabled by '!')
            line = self._prefix + self._key
        else:
            # integer, double, or string parameter
            line = self._prefix + self._key + Keyword.fourspaces + str(self._value)

        linelength = len(line)
        return linelength, line


class BooleanKeyword(Keyword):
    """Chemkin boolean keyword."""

    def __init__(self, phrase: str):
        """Set up a Chemkin keyword with a boolean parameter."""
        value = True
        super().__init__(phrase, value, "bool")


class IntegerKeyword(Keyword):
    """A Chemkin integer keyword."""

    def __init__(self, phrase: str, value: int = 0):
        """Set up a Chemkin keyword with an integer parameter."""
        super().__init__(phrase, value, "int")


class RealKeyword(Keyword):
    """A Chemkin real keyword."""

    def __init__(self, phrase: str, value: float = 0.0e0):
        """Set up a Chemkin keyword with a real number parameter."""
        super().__init__(phrase, value, "float")


class StringKeyword(Keyword):
    """A Chemkin string keyword."""

    def __init__(self, phrase: str, value: str = ""):
        """Set up a Chemkin keyword with a string parameter."""
        if len(value) <= 0:
            msg = [Color.PURPLE, "no string parameter given", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        super().__init__(phrase, value, "str")
