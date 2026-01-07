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

r""".. _ref_multiple_mechanism:

=============================
Work with multiple mechanisms
=============================

PyChemkin can facilitate multiple mechanisms in one project. However, only one
chemistry set can be active at a time. This example shows how to use
the ``activate()`` method to switch between multiple chemistry sets (mechanisms)
in the same Python project. You can use this method to compare results from
two different mechanisms, such as the *base* and *reduced* mechanisms.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_preprocess.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)


################################
# Create the first chemistry set
# ==============================
# The first mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir

# specify the mechanism input files
# including the full file path is recommended
chemfile = str(mechanism_dir / "grimech30_chem.inp")
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
tranfile = str(mechanism_dir / "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
My1stMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")


####################################
# Preprocess the first chemistry set
# ==================================

# preprocess the mechanism files
ierror = My1stMech.preprocess()
print()
if ierror != 0:
    # encountered error during preprocessing
    print(f"Preprocessing error encountered. Code = {ierror:d}.")
    print(f"See the summary file {My1stMech.summaryfile} for details.")
    exit()
else:
    # Display the basic mechanism information
    print(Color.GREEN + "Preprocessing succeeded.", end=Color.END)
    print("Mechanism information:")
    print(f"Number of elements = {My1stMech.mm:d}.")
    print(f"Number of gas species = {My1stMech.kk:d}.")
    print(f"Number of gas reactions = {My1stMech.ii_gas:d}.")


######################
# Set up a gas mixture
# ====================
# Set up a gas mixture based on the species in this chemistry set.
# Create a gas mixture named ``mymixture1`` based on the ``My1stMech`` chemistry set.
# The species mole fractions/ratios are given in the ``recipe`` format. The ``X``
# property is used because the mole fractions are given.

mymixture1 = ck.Mixture(My1stMech)
# set mixture temperature [K]
mymixture1.temperature = 1000.0
# set mixture pressure [dynes/cm2]
mymixture1.pressure = ck.P_ATM
# use the "X" property to specify the molar compositions of the mixture
mymixture1.x = [("CH4", 0.1), ("O2", 0.21), ("N2", 0.79)]


####################################
# Perform an equilibrium calculation
# ==================================
# The equilibrium state is stored as a mixture named ``equil_mix1_hp``. You can get
# the equilibrium temperature from the ``temperature`` property of this mixture.
# You can learn more about the PyChemkin ``equilibrium()`` method by either typing
# ``ck.help("equilibrium")`` at the Python prompt or uncommenting the following line.

# ck.help("equilibrium")

# find the constrained H-P equilibrium state of ``mymixture1``
equil_mix1_hp = ck.equilibrium(mymixture1, opt=5)
# print the equilibrium temperature
print(f"Equilibrium temperature of mymixture1: {equil_mix1_hp.temperature} [K]")


#################################
# Create the second chemistry set
# ===============================
# The second mechanism is the ``C2 NOx`` mechanism in the ``/reaction/data``
# directory. Here, a different process for setting up a chemistry set is used.
# The chemistry set instance is created before the mechanism files are specified.
# You can make changes to the files to include in the chemistry set before running
# the preprocessing step.
#
# The ``C2 NOx`` mechanism file, in addition to the reactions, contains the
# thermodynamic and transport data of all species in the mechanism. Thus,
# you must only specify the mechanism file, that is, ``chemfile``. If your simulation
# requires transport properties, you must use the ``preprocess_transportdata()``
# method to tell the preprocessor to also include the transport data.

# set the second mechanism directory (the default Chemkin mechanism data directory)
mechanism_dir = data_dir
# create a chemistry set based on "C2_NOx" using an alternative method
My2ndMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all necessary thermodynamic and transport data
# thus, there is no need to specify thermodynamic and transport data files
My2ndMech.chemfile = str(mechanism_dir / "C2_NOx_SRK.inp")

# direct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
My2ndMech.preprocess_transportdata()


#####################################
# Preprocess the second chemistry set
# ===================================
# The ``C2 NOx`` mechanism also includes information about the *Soave* cubic
# Equation of State (EOS) for real-gas applications. The PyChemkin preprocessor
# indicates the availability of the real-gas model in the chemistry set processed.
# For example, during preprocessing, this is printed:
# ``real-gas cubic EOS 'Soave' is available``. As soon as the second chemistry set
# is preprocessed successfully, it becomes the active chemistry set of the project.
# The first chemistry set, ``My1stMech``, is pushed to the background.

# preprocess the second mechanism files
ierror = My2ndMech.preprocess()
print()
if ierror != 0:
    # encountered error during preprocessing
    print(f"Preprocessing error encountered. Code = {ierror:d}.")
    print(f"See the summary file {My2ndMech.summaryfile} for details.")
    exit()
else:
    # Display the basic mechanism information
    print(Color.GREEN + "Preprocessing succeeded.", end=Color.END)
    print("Mechanism information:")
    print(f"Number of elements = {My2ndMech.mm:d}.")
    print(f"Number of gas species = {My2ndMech.kk:d}.")
    print(f"Number of gas reactions = {My2ndMech.ii_gas:d}.")


################################################################################
# Set up the second gas mixture based on the species in the second chemistry set
# ==============================================================================
# Create a gas mixture named ``mymixture2`` based on the ``My2ndMech`` chemistry set.
mymixture2 = ck.Mixture(My2ndMech)
# set mixture temperature [K]
mymixture2.temperature = 500.0
# set mixture pressure [dynes/cm2]
mymixture2.pressure = 2.0 * ck.P_ATM
# set mixture molar composition
mymixture2.x = [("H2", 0.02), ("O2", 0.2), ("N2", 0.8)]


###############################
# Run the detonation calculation
# ==============================
# You can now compute the detonation wave speed with the ``mymixture2`` gas mixture.
# The ``CJ_mix2`` represents the mixture at the Chapman-Jouguet state. The
# speed of sound and the detonation wave speed are returned in
# the ``speeds_mix2`` tuple.

speeds_mix2, cj_mix2 = ck.detonation(mymixture2)
#  print the detonation calculation results
print(f"Detonation mymixture2 temperature: {cj_mix2.temperature} [K]")
print(f"Detonation wave speed = {speeds_mix2[1] / 100.0} [m/sec]")


###################################################
# Switch to the first chemistry set and gas mixture
# =================================================
# Use the ``activate()`` method to reactivate the ``My1stMech`` chemistry set and
# the ``mymixture1`` gas mixture.

My1stMech.activate()


################################
# Run the detonation calculation
# ==============================
# You can now compute the detonation wave speed with the ``mymixture1`` gas mixture.
# The ``CJ_mix1`` represents the mixture at the Chapman-Jouguet state. The
# speed of sound and the detonation wave speed are returned in
# the ``speeds_mix1`` tuple.
#
# .. note::
#   The ``mymixture1`` and ``mymixture2`` gas mixtures have different
#   initial conditions.

speeds_mix1, cj_mix1 = ck.detonation(mymixture1)
# print the detonation calculation results
print(f"Detonation 'mymixture1' temperature: {cj_mix1.temperature} [K]")
print(f"Detonation wave speed = {speeds_mix1[1] / 100.0} [m/sec]")
