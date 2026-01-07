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

"""Chemkin module for Python core."""

from ctypes import c_int
import inspect
from pathlib import Path
import platform

# import Chemkin-CFD-API
# import all commonly used constants and methods
# so the users can have easy access to these resources
from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.chemistry import (
    Chemistry as Chemistry,
    chemkin_version as chemkin_version,
    done as done,
    set_verbose as set_verbose,
    verbose as verbose,
)
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import (
    AVOGADRO as AVOGADRO,
    BOLTZMANN as BOLTZMANN,
    ERGS_PER_CALORIE as ERGS_PER_CALORIE,
    ERGS_PER_JOULE as ERGS_PER_JOULE,
    JOULES_PER_CALORIE as JOULES_PER_CALORIE,
    P_ATM as P_ATM,
    P_TORRS as P_TORRS,
    R_GAS as R_GAS,
    R_GAS_CAL as R_GAS_CAL,
    Air as Air,
    water_heat_vaporization as water_heat_vaporization,
)
from ansys.chemkin.core.info import (
    help as help,
    keyword_hints as keyword_hints,
    manuals as manuals,
    phrase_hints as phrase_hints,
    setup_hints,
    show_equilibrium_options as show_equilibrium_options,
    show_ignition_definitions as show_ignition_definitions,
    show_realgas_usage as show_realgas_usage,
)
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import (
    Mixture as Mixture,
    adiabatic_mixing as adiabatic_mixing,
    cal_mixture_temperature_from_enthalpy as cal_mixture_temperature_from_enthalpy,
    calculate_equilibrium as calculate_equilibrium,
    detonation as detonation,
    equilibrium as equilibrium,
    interpolate_mixtures as interpolate_mixtures,
    isothermal_mixing as isothermal_mixing,
)
from ansys.chemkin.core.realgaseos import (
    check_realgas_status as check_realgas_status,
    set_current_pressure as set_current_pressure,
)

# show ansys (chemkin) version number
msg = [
    Color.YELLOW,
    "Chemkin version number =",
    str(chemkin_version()),
    Color.END,
]
this_msg = Color.SPACE.join(msg)
logger.info(this_msg)
# get ansys installation location
ansys_dir = ck_wrapper._ansys_dir
ansys_version = ck_wrapper._ansys_ver

if platform.system() == "Windows":
    _chemkin_platform = "win64"
else:
    _chemkin_platform = "linuxx8664"

# get chemkin installation location
ck_name = "chemkin." + _chemkin_platform
_chemkin_root = Path(ansys_dir) / "reaction" / ck_name
chemkin_dir = _chemkin_root.resolve()

# set default units to cgs
unit_code = c_int(1)
ierror = ck_wrapper.chemkin.KINSetUnitSystem(unit_code)

# chemkin module home directory
frm = inspect.currentframe()
if frm is not None:
    _chemkin_module_path = Path(inspect.getfile(frm))
    pychemkin_dir = _chemkin_module_path.parent
# set up Chemkin keyword help data
setup_hints()
