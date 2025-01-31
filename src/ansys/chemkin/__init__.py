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

"""
Chemkin module for Python
core
"""

from ctypes import c_int
import inspect
import os
from pathlib import Path
import platform

# import Chemkin-CFD-API
# import all commonly used constants and methods
# so the users can have easy access to these resources
from . import chemkin_wrapper as ck_wrapper
from .chemistry import Chemistry, done, set_verbose, verbose
from .color import Color
from .constants import (
    Patm,
    RGas,
    RGas_Cal,
    avogadro,
    boltzmann,
    ergs_per_calorie,
    ergs_per_joule,
    joules_per_calorie,
)
from .info import (
    help,
    keyword_hints,
    manuals,
    phrase_hints,
    setup_hints,
    show_equilibrium_options,
    show_ignition_definitions,
    show_realgas_usage,
)
from .logger import logger
from .mixture import (
    Mixture,
    adiabatic_mixing,
    calculate_equilibrium,
    calculate_mixture_temperature_from_enthalpy,
    detonation,
    equilibrium,
    interpolate_mixtures,
    isothermal_mixing,
)
from .realgaseos import check_realgas_status, set_current_pressure

# show ansys (chemkin) version number
msg = [
    Color.YELLOW,
    "Chemkin version number =",
    str(ck_wrapper._ansys_ver),
    Color.END,
]
this_msg = Color.SPACE.join(msg)
logger.info(this_msg)
# get ansys installation location
ansys_dir = str(ck_wrapper._ansys_dir)
ansys_version = ck_wrapper._ansys_ver

_chemkin_platform = None
if platform.system() == "Windows":
    _chemkin_platform = "win64"
else:
    _chemkin_platform = "linuxx8664"

# get chemkin installation location
_chemkin_root = os.path.join(ansys_dir, "reaction", "chemkin." + _chemkin_platform)
chemkin_dir = Path(_chemkin_root)

# set default units to cgs
unit_code = c_int(1)
iError = ck_wrapper.chemkin.KINSetUnitSystem(unit_code)

# chemkin module home directory
frm = inspect.currentframe()
if frm is not None:
    _chemkin_module_path = os.path.dirname(inspect.getfile(frm))
    pychemkin_dir = Path(_chemkin_module_path)
# set up Chemkin keyword help data
setup_hints()
