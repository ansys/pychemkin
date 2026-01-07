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

"""Color and print format utilities."""


class Color:
    """Define colors used by PyChemkin for printing text messages."""

    # color codes
    RED = "\033[91m"
    MAGENTA = "\033[035m"
    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    BLUE = "\033[94m"
    YELLOW = "\033[93m"
    WHITE = "\033[037m"
    BOLD = "\033[1m"
    ULINE = "\033[4m"
    END = "\n \033[0m"
    SPACE = " "
    SPACEx6 = "      "
    # message colors
    msg_modes = {
        "normal": WHITE,
        "info": YELLOW,
        "warning": MAGENTA,
        "error": PURPLE,
        "critical": RED,
        "ok": GREEN,
    }
    # logger levels
    log_modes = {
        "normal": 0,
        "info": 20,
        "warning": 30,
        "error": 40,
        "critical": 50,
        "ok": 10,
    }

    @staticmethod
    def ckprint(mode: str, msg: list = []):
        """Customize text messages."""
        """
        Customize text messages.

        Parameters
        ----------
            mode: string, {"normal, "info", "warning, "error", "fatal", "ok"},
            default = ""
                message mode/type
            msg: list of strings
                the message to be printed

        """
        # check message to be printed
        if len(msg) <= 0:
            # no message
            return
        # set color
        color = Color.msg_modes.get(mode.lower(), "WHITE")
        # compile the message
        message = Color.SPACE.join(msg)
        # print
        print(color + f"** {message}", end=Color.END)
