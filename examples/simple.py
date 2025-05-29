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

==========================
Set up a PyChemkin project
==========================

This example shows how to set up a PyChemkin project and start to use PyChemkin features.
You begin by importing the PyChemkin package, which is named ``ansys.chemkin``, and optionally
the PyChemkin logger package, which is named ``ansys.chemkin.logger``.

Next, you use the ``Chemistry()`` method to instantiate and preprocess a chemistry set.
This requires specifying the file names (with full file paths) of the mechanism file,
thermodynamic data file, and optional transport data file.

Once the chemistry set instance is instantiated, you use the ``preprocess()`` method to
preprocess the mechanism data. You can then start to use PyChemkin features to build
your simulation.

.. note::
   Running PyChemkin requires Ansys 2025 R2 or later.
"""

###############################################
# Import PyChemkin packages
# =========================
# Import the ``ansys.chemkin`` package to start using PyChemkin in the project.
# Importing the ``ansys.chemkin.logger`` package is optional.

import os

import ansys.chemkin  # import PyChemkin
from ansys.chemkin.logger import logger

#######################################################
# Set up the file paths of the mechanism and data files
# =====================================================
# The ``ansys_dir`` variable on your local computer specifies the location of your Ansys installation.
#
# Use ``os.path`` methods to construct the file names with the full paths.
#
# Assign the files to the corresponding chemistry set arguments:
#
# - ``mech_file``: Mechanism input file
# - ``therm_file``: Thermodynamic data file
# - ``tran_file``: Transport data file (optional)
#

# create GRI 3.0 mechanism from the data directory
mechanism_dir = os.path.join(ansys.chemkin.ansys_dir, "reaction", "data")
# set up mechanism file names
mech_file = os.path.join(mechanism_dir, "grimech30_chem.inp")
therm_file = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tran_file = os.path.join(mechanism_dir, "grimech30_transport.dat")

###############################
# Instantiate the chemistry set
# =============================
# Use the ``Chemistry()`` method to instantiate a chemistry set named ``GasMech``.

GasMech = ansys.chemkin.Chemistry(
    chem=mech_file, therm=therm_file, tran=tran_file, label="GRI 3.0"
)

##############################
# Preprocess the chemistry set
# ============================
# Preprocess the GRI 3.0 chemistry set.

status = GasMech.preprocess()

#################################
# Check the preprocessing status
# ===============================

if status != 0:
    # failed
    print(f"Preprocessing: Error encountered. Code = {status:d}.")
    print(f"See the summary file {GasMech.summaryfile} for details.")
    logger.error("Preprocessing failed.")
    exit()

#################################
# Start to use PyChemkin features
# ===============================
# Start to use PyChemkin features to build your simulation.
# For example, using the ``Mixture()`` method, create an ``air`` mixture based
# on the ``GasMech`` chemistry set.

# create 'air' mixture based on 'GasMech' chemistry set
air = ansys.chemkin.Mixture(GasMech)
# set 'air' condition
# mixture pressure in [dynes/cm2]
air.pressure = 1.0 * ansys.chemkin.Patm
# mixture temperature in [K]
air.temperature = 300.0
# mixture composition in mole fractions
air.X = [("O2", 0.21), ("N2", 0.79)]

##########################################################
# Print the properties of the air mixture for verification
# ========================================================
#
# .. note::
#
#   - The default units of temperature and pressure are [K] and [dynes/cm\ :sup:`2`\ ], respectively.
#   - The constant ``Patm`` is a conversion multiplier for pressure.
#   - Transport property methods such as ``mixture_viscosity()`` require transport data. You
#     must include the transport data file when creating the chemistry set.
#

# print pressure and temperature of the `air` mixture
print(f"Pressure    = {air.pressure/ansys.chemkin.Patm} [atm]")
print(f"Temperature = {air.temperature} [K]")
# print the 'air' composition in mass fractions
air.list_composition(mode="mass")
# get 'air' mixture density [g/cm3]
print(f"Mixture density   = {air.RHO} [g/cm3]")
# get 'air' mixture viscosity [g/cm-sec] or [poise]
print(f"Mixture viscosity = {air.mixture_viscosity()*100.0} [cP]")
