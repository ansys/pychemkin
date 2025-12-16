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
.. _ref_load_mechanism:

================================
Preprocess a gas-phase mechanism
================================

Before you can run a Chemkin simulation, you must load and preprocess the reaction
mechanism data. The mechanism data consists of a reaction mechanism input file,
a thermodynamic data file, and an optional transport data file.

This example shows how to instantiate and preprocess a chemistry set, which is
always the first task in running any PyChemkin simulation.

PyChemkin allows several chemistry sets to coexist in the same Python project.
However, only one chemistry set can be active at a time.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_preprocess.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os

# import PyChemkin packages
import ansys.chemkin as ck
from ansys.chemkin import Color
from ansys.chemkin.logger import logger

# check the working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)

# set PyChemkin verbose mode
ck.set_verbose(True)


########################
# Create a chemistry set
# ======================
# The first mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir

# set mechanism input files
# including the full file path is recommended
# the gas-phase reaction mechanism file (GRI 3.0)
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
# the thermodynamic data file
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# the transport data file
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")

# create a chemistry set based on the GRI 3.0 methane combustion mechanism
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")


##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

# display the preprocess status
print()
if iError != 0:
    # When a non-zero value is returned from the process, check the text output files,
    # "chem.out," "tran.out," and "summary.out," for potential error messages about the mechanism data.
    print(f"Preprocessing error encountered. Code = {iError:d}.")
    print(f"See the summary file {MyGasMech.summaryfile} for details.")
    exit()
else:
    Color.ckprint("OK", ["Preprocessing succeeded.", "!!!"])
    print("Mechanism information:")
    print(f"Number of elements = {MyGasMech.number_elements:d}.")
    print(f"Number of gas species = {MyGasMech.number_species:d}.")
    print(f"Number of gas reactions = {MyGasMech.number_gas_reactions:d}.")


#########################################
# Display the basic mechanism information
# =======================================

print(f"\nElement and species information of mechanism {MyGasMech.label}")
print("=" * 50.0)

# extract element symbols as a list
elelist = MyGasMech.element_symbols
# get atomic masses as numpy 1D double array
AWT = MyGasMech.atomic_weight
# print element information
for k in range(len(elelist)):
    print(f"Element number {k+1:3d}: {elelist[k]:16}. Mass = {AWT[k]:f}.")

print("=" * 50)

# extract gas species symbols as a list
specieslist = MyGasMech.species_symbols
# get species molecular masses as numpy 1D double array
WT = MyGasMech.species_molar_weight
# print gas species information
for k in range(len(specieslist)):
    print(f"Species number {k+1:3d}: {specieslist[k]:16}. Mass = {WT[k]:f}.")
print("=" * 50)


###############################
# Create a second chemistry set
# =============================
# The second mechanism to load into the project is the C2-NOx mechanism for
# the combustion of C1-C2 hydrocarbons. This mechanism differs from the GRI mechanism
# in the sense that it is self-contained, that is, the thermodynamic and transport data
# of all species is included in the ``C2_NOx_SRK.inp`` mechanism input file, which sits in
# the same reaction data directory as the GRI mechanism input files.
#
# Use the same steps that were used to instantiate the first chemistry set to
# process any set of reaction mechanism files. Here, use a slightly different procedure
# to instantiate the second chemistry set. After you have created it, specify
# the reaction mechanism files one by one. The reaction mechanism
# file in this case contains all the necessary thermodynamic and transport data. Thus, there
# is no need to specify thermodynamic and transport data files. However, an additional step is
# required to instruct the preprocessor to include the transport data.

# set the second mechanism directory (the default Chemkin mechanism data directory)
mechanism_dir = data_dir

# create a chemistry set based on C2_NOx using an alternative method
My2ndMech = ck.Chemistry(label="C2 NOx")

# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# thus, there is no need to specify thermodynamic and transport data files
My2ndMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

# instruct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
My2ndMech.preprocess_transportdata()


#####################################
# Preprocess the second chemistry set
# ===================================

# preprocess the second mechanism file
iError = My2ndMech.preprocess()

# display the preprocess status
print()
if iError != 0:
    # When a non-zero value is returned from the process, check the text output files,
    # "chem.out," "tran.out," and "summary.out," for potential error messages about the mechanism data.
    print(f"Preprocessing error encountered. Code = {iError:d}.")
    print(f"See the summary file {My2ndMech.summaryfile} for details.")
    exit()
else:
    Color.ckprint("OK", ["Preprocessing succeeded.", "!!!"])
    print("Mechanism information:")
    print(f"Number of elements = {My2ndMech.MM:d}.")
    print(f"Number of gas species = {My2ndMech.KK:d}.")
    print(f"Number of gas reactions = {My2ndMech.IIGas:d}.")


#########################################
# Display the basic mechanism information
# =======================================

print(f"\nElement and species information of mechanism {My2ndMech.label}")
print("=" * 50.0)

# extract element symbols as a list
elelist = My2ndMech.element_symbols
# get atomic masses as numpy 1D double array
AWT = My2ndMech.AWT
# print element information
for k in range(len(elelist)):
    print(f"Element # {k+1:3d}: {elelist[k]:16}. Mass = {AWT[k]:f}.")

print("=" * 50)

# extract gas species symbols as a list
specieslist = My2ndMech.species_symbols
# get species molecular masses as numpy 1D double array
WT = My2ndMech.WT
# print gas species information
for k in range(len(specieslist)):
    print(f"Species # {k+1:3d}: {specieslist[k]:16}. Mass = {WT[k]:f}.")

print("=" * 50)
