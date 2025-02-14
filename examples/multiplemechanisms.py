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
.. _ref_multiple_mechanism:

=======================================================
Working with multiple mechanisms in a PyChemkin project
=======================================================

PyChemkin can facilitate multiple mechanisms in the one project; however, only one active
``Chemistry Set`` object is allowed at a time. This tutorial demonstrates
how to switch between multiple ``Chemistry Set`` objects (mechanisms) in the same python project by
using the ``activate`` method. You can use this feature to compare the results from two different mechanisms
such as the *master* and the *reduced* mechanisms.
"""

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.logger import logger

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)


####################################################
# Create an instance of the first  ``Chemistry Set``
# ==================================================
# The first mechanism loaded is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir

# specify the mechanism input files
# inclusion of the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
My1stMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")


###########################################
# Pre-process the GRI 3.0 ``Chemistry Set``
# =========================================

# preprocess the mechanism files
iError = My1stMech.preprocess()
print()
if iError != 0:
    # encountered error during preprocessing
    print(f"PreProcess: error encountered...code = {iError:d}")
    print(f"see the summary file {My1stMech.summaryfile} for details")
    exit()
else:
    # Display the basic mechanism information
    print(Color.GREEN + "PreProcess success!!", end=Color.END)
    print("mechanism information:")
    print(f"number of elements = {My1stMech.MM:d}")
    print(f"number of gas species = {My1stMech.KK:d}")
    print(f"number of gas reactions = {My1stMech.IIGas:d}")


#############################################################################
# Set up a gas mixture based on the species in this GRI 3.0 ``Chemistry Set``
# ===========================================================================
# Create gas mixture instance ``mymixture1`` based on ``My1stMech`` ``Chemistry Set``.
# The species mole fractions/ratios are given in the  ``recipe`` format, and the ``X``
# method is used here because of the mole fractions are given.

mymixture1 = ck.Mixture(My1stMech)
# set mixture temperature [K]
mymixture1.temperature = 1000.0
# set mixture pressure [dynes/cm2]
mymixture1.pressure = ck.Patm
# use the "X" property to specify the molar compositions of the mixture
mymixture1.X = [("CH4", 0.1), ("O2", 0.21), ("N2", 0.79)]


########################################
# Perform an ``equilibrium`` calculation
# ======================================
# The equilibrium state will be stored as a ``Mixure``: ``equil_mix1_HP``. You can get
# the equilibrium temperature from the ``temperature`` property of ``equil_mix1_HP``.
# You can learn more about the PyChemkin ``equilibrium`` method by typing
# *ck.help("equilibrium")* at the python prompt, or uncomment the line below.

# ck.help("equilibrium")

# Find the *constrained H-P* equilibrium state of ``mymixture1``.
equil_mix1_HP = ck.equilibrium(mymixture1, opt=5)
# Print the equilibrium temperature
print(f"equilibrium temperature of mymixture1 : {equil_mix1_HP.temperature} [K]")


####################################################
# Create an instance of the second ``Chemistry Set``
# ==================================================
# The second mechanism is the 'C2 NOx' mechanism from the default *"/reaction/data"*
# directory. A different ``Chemistry Set`` setup process is employed here.
# The ``Chemistry Set`` instance is created before the mechanism files are specified.
# You can make changes to the files to be included in the ``Chemistry Set`` before running
# the ``preprocess`` step.
# The 'C2 NOx' mechanism file, in addition to the reactions, also contains the thermodynamic
# and the transport data of all species in the mechanism. In this case, you only need to specify
# the mechanism file, that is, ``chemfile``. If your simulation requires the transport properties, you
# must use the ``preprocess_transportdata`` method to tell the preprocessor to also include the transport data.

# set the 2nd mechanism directory (the default chemkin mechanism data directory)
mechanism_dir = data_dir
# create a chemistry set based on C2_NOx using an alternative method
My2ndMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
My2ndMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

# direct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
My2ndMech.preprocess_transportdata()


##########################################
# Pre-process the C2 NOx ``Chemistry Set``
# ========================================
# The 'C2 NOx' mechanism also includes information about the *Soave* cubic
# Equation of State (EOS) for the real-gas applications. PyChemkin preprocessor
# will indicate the availability of the real-gas model in the ``Chemistry Set`` processed.
# For example, you will see the print-out *"real-gas cubic EOS 'Soave' is available"* during
# the preprocess.
# As soon as the second ``Chemistry Set`` is preprocessed successfully, it becomes the **active**
# ``Chemistry Set`` of the project. The first ``Chemistry Set``, ``My1stMech`` is pushed to the background.

# preprocess the 2nd mechanism files
iError = My2ndMech.preprocess()
print()
if iError != 0:
    # encountered error during preprocessing
    print(f"PreProcess: error encountered...code = {iError:d}")
    print(f"see the summary file {My2ndMech.summaryfile} for details")
    exit()
else:
    # Display the basic mechanism information
    print(Color.GREEN + "PreProcess success!!", end=Color.END)
    print("mechanism information:")
    print(f"number of elements = {My2ndMech.MM:d}")
    print(f"number of gas species = {My2ndMech.KK:d}")
    print(f"number of gas reactions = {My2ndMech.IIGas:d}")


#####################################################################################
# Set up the second gas mixture based on the species in the C2 NOx ``Chemistry Set``
# ===================================================================================
# Create gas mixture instance ``mymixture2`` based on ``My2ndMech`` ``Chemistry Set``.
mymixture2 = ck.Mixture(My2ndMech)
# set mixture temperature [K]
mymixture2.temperature = 500.0
# set mixture pressure [dynes/cm2]
mymixture2.pressure = 2.0 * ck.Patm
# set mixture molar composition
mymixture2.X = [("H2", 0.02), ("O2", 0.2), ("N2", 0.8)]


########################################
# Perform a ``detonation`` calculation
# ======================================
# Now you can compute detonation wave speed with ``mymixture2``. ``CJ_mix2`` represents the mixture at
# the Chapman-Jouguet state, and the *speed of sound* and the *detonation wave speed* are returned in
# the tuple ``speeds_mix2``.

speeds_mix2, CJ_mix2 = ck.detonation(mymixture2)
#  print the detonation calculation results
print(f"detonation mymixture2 temperature: {CJ_mix2.temperature} [K]")
print(f"detonation wave speed = {speeds_mix2[1]/100.0} [m/sec]")


#############################################################################################
# Switch back to the first ``Chemistry Set`` ``My1stMech`` and the gas mixture ``mymixture1``
# ===========================================================================================
# Use the ``activate`` method to re-activate ``My1stMech`` and ``mymixture1``.
My1stMech.activate()


########################################
# Perform a ``detonation`` calculation
# ======================================
# Now you can compute detonation wave speed with ``mymixture1``. ``CJ_mix1`` represents the mixture at
# the Chapman-Jouguet state, and the *speed of sound* and the *detonation wave speed* are returned in
# the tuple ``speeds_mix1``.
#
# .. note::
#   ``mymixture1`` and ``mymixture2`` have different initial conditions.

speeds_mix1, CJ_mix1 = ck.detonation(mymixture1)
#  print the detonation calculation results
print(f"detonation mymixture1 temperature: {CJ_mix1.temperature} [K]")
print(f"detonation wave speed = {speeds_mix1[1]/100.0} [m/sec]")
