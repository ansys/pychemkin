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
.. _ref_simple:

=====================================================
Tasks required at the start of all PyChemkin projects
=====================================================

This tutorial walks you through the start-up procedures that are necessary to every
PyChemkin project.

Firstly, the PyChemkin package must be imported. The official package name of PyChemkin is ``ansys.chemkin``.
You can also import the PyChemkin logger ``ansys.chemkin.logger``, but this step is optional.
Before performing any calculations, you must also **instantiate** and **preprocess**
a ``Chemistry Set``. The ``Chemistry Set`` instance is created by using the ``Chemistry`` method to
which you specify the *file names* (with full file path) of the *mechanism file*, the
*thermodynamic data file*, and the optional *transport data file*. Once the ``Chemistry Set`` is instantiated,
use the ``preprocess`` method to preprocess the mechanism data. Once these steps are finished successfully,
you may start to build your PyChemkin simulation project.
"""

###############################################
# Import PyChemkin package and start the logger
# =============================================
# Import the ``ansys.chemkin`` package to start using PyChemkin in the project.
# The ``ansys.chemkin.logger`` is optional.

import os

import ansys.chemkin  # import PyChemkin
from ansys.chemkin.logger import logger

#################################################################
# Set up the file paths of the mechanism input and the data files
# ===============================================================
# ``ansys_dir`` is the **Ansys** installation on your local computer.
#
# Use ``os.path`` methods to construct the file names with full path.
#
# Assign the files to the corresponding ``Chemistry Set`` arguments
#
#   | ``chem`` : mechanism input file
#   | ``therm``: thermodynamic data file
#   | ``tran`` : transport data file (optional)
#
#
# .. note::
#   The minimum version to run PyChemkin is **Ansys 2025 Release 2**.
#

# create GRI 3.0 mechanism from the data directory.
mechanism_dir = os.path.join(ansys.chemkin.ansys_dir, "reaction", "data")
# set up mechanism file names
mech_file = os.path.join(mechanism_dir, "grimech30_chem.inp")
therm_file = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tran_file = os.path.join(mechanism_dir, "grimech30_transport.dat")

###################################
# Instantiate the ``Chemistry Set``
# =================================
# Use the ``Chemistry`` method to instantiate ``Chemistry Set`` ``GasMech``.

GasMech = ansys.chemkin.Chemistry(
    chem=mech_file, therm=therm_file, tran=tran_file, label="GRI 3.0"
)

###########################################
# Pre-process the GRI 3.0 ``Chemistry Set``
# =========================================
# Preprocess the Chemistry Set.
status = GasMech.preprocess()

#################################
# Verify the preprocessing status
# ===============================
# Check preprocess status.
if status != 0:
    # failed
    print(f"PreProcess: error encountered...code = {status:d}")
    print(f"see the summary file {GasMech.summaryfile} for details")
    logger.error("PreProcess failed")
    exit()

#################################
# Start to use PyChemkin features
# ===============================
# For example, you can create an "air" mixture based on ``GasMech``
# by using the ``Mixutre`` method.

# Create Mixture 'air' based on 'GasMech'
air = ansys.chemkin.Mixture(GasMech)
# set 'air' condition
# mixture pressure in [dynes/cm2]
air.pressure = 1.0 * ansys.chemkin.Patm
# mixture temperature in [K]
air.temperature = 300.0
# mixture composition in mole fractions
air.X = [("O2", 0.21), ("N2", 0.79)]

##############################################################
# Print the properties of the ``air`` mixture for verification
# ============================================================
#
# .. note::
#   The default units of temperature and pressure are [K] and [dynes/cm\ :sup:`2`\ ], respectively.
#
#   The constant ``Patm`` is a conversion multiplier for pressure.
#
#
# .. note::
#   Transport property methods such as ``mixture_viscosity`` require *transport data*; you
#   must include the ``tran`` data file when creating the ``Chemistry Set``.
#

# print pressure and temperature of the `air` mixture
print(f"pressure    = {air.pressure/ansys.chemkin.Patm} [atm]")
print(f"temperature = {air.temperature} [K]")
# print the 'air' composition in mass fractions
air.list_composition(mode="mass")
# get 'air' mixture density [g/cm3]
print(f"the mixture density   = {air.RHO} [g/cm3]")
# get 'air' mixture viscosity [g/cm-sec] or [poise]
print(f"the mixture viscosity = {air.mixture_viscosity()*100.0} [cP]")
