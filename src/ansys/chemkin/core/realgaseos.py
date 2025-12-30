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

"""Real-gas cubic EOS model."""

from ctypes import c_double, c_int

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper


def check_realgas_status(chem_index: int) -> bool:
    """Check whether the real-gas cubic EOS is active."""
    """
    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    Returns
    -------
        status: boolean
            the activation status of the Chemkin real-gas model

    """
    # initialization assuming the real-gas EOS is not ON
    status = False

    chemset_index = c_int(chem_index)
    mode = c_int(0)
    ierr = ck_wrapper.chemkin.KINRealGas_CheckRealGasStatus(chemset_index, mode)
    if ierr == 0:
        status = mode.value == 1
    return status


def set_current_pressure(chem_index: int, pressure: float) -> int:
    """Set gas mixture pressure for real-gas EOS calculations."""
    """
    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set
        pressure: double
            gas pressure [dynes/cm2]

    Returns
    -------
        ierr: integer
            error code

    """
    # convert variables
    chemset_index = c_int(chem_index)
    p = c_double(pressure)
    ierr = ck_wrapper.chemkin.KINRealGas_SetCurrentPressure(chemset_index, p)
    return ierr
