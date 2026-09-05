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

r""".. _ref_species_properties:

===============================
Evaluate gas species properties
===============================

PyChemkin inherits many Ansys Chemkin utilities for evaluating
species properties in a mechanism. Most tools for calculating
the species properties are part of the PyChemkin
``Chemistry()`` methods. This example shows how use some of
these methods to evaluate species thermodynamic and
transport properties.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_species_properties.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True


########################
# Create a chemistry set
# ======================
# The mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir

# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusding the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
MyGasMech.tranfile = str(mechanism_dir / "grimech30_transport.dat")


##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()


#########################################
# Display the basic mechanism information
# =======================================
# Display the element and species information from the GRI 3.0 mechanism.
# You can get the entire list of the elements and the gas species in
# the mechanism by using the ``element_symbols`` and
# ``species_symbols`` lists, respectively.
#
# Use the ``get_specindex()`` method to get the species index of a species
# in the mechanism. Use the ``SpeciesComposition()`` method to check
# the elemental composition of a gas species.
#
# .. note::
#   Both the species name/symbol and element name/symbol are case-sensitive.

# extract element symbols as a list
elelist = MyGasMech.element_symbols
# extract gas species symbols as a list
specieslist = MyGasMech.species_symbols
# list of gas species of interest
plotspeclist = ["CH4", "O2", "N2"]
# find elemental compositions of selected species
print(" ")
for s in plotspeclist:
    species_index = MyGasMech.get_specindex(s)
    print("species " + specieslist[species_index])
    print("elemental composition")
    for elem_index in range(MyGasMech.mm):
        num_elem = MyGasMech.species_composition(elem_index, species_index)
        print(f"    {elelist[elem_index]:>4}: {num_elem:2d}")
    print("=" * 10)
print()


##################################
# Plot selected species properties
# ================================
# Calculate and plot the properties of CH\ :sub:`4`\ , O\ :sub:`2`\ , and N\ :sub:`2`
# against the temperature. Use these property methods:
#
# - ``species_cv``: Species specific heat capacity at constant volume [erg/mol-K]
# - ``SpeciesCond``: Species thermal conductivity [erg/cm-K-sec]
# - ``SpeciesDiffusionCoeffs``: Binary diffusion coefficients between species pairs
#   [cm\ :sup:`2`\ /sec]
#

# plot Cv and thermal conductivity values at different temperatures for
# selected gas species
#
plt.figure(figsize=(12, 6))
# temperature increment
dtemp = 20.0
# number of property data points
points = 100
# curve attributes
curvelist = ["g", "b--", "r:"]


##########################################
# Prepare the species Cv data for plotting
# ========================================

# create arrays
# specify specific heat capacity at constant volume data
cv = np.zeros(points, dtype=np.double)
# temperature data
t = np.zeros(points, dtype=np.double)
# start plotting loop #1
k = 0
# loop over selected gas species
for s in plotspeclist:
    # starting temperature at 300 [K]
    temp = 300.0
    # loop over temperature data points
    for i in range(points):
        heatcapacity = MyGasMech.species_cv(temp)
        index = MyGasMech.get_specindex(s)
        t[i] = temp
        # convert ergs to joules
        cv[i] = heatcapacity[index] / ck.ERGS_PER_JOULE
        temp += dtemp
    plt.subplot(121)
    plt.plot(t, cv, curvelist[k])
    k += 1
# plot Cv versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Cv [J/mol-K]")
plt.legend(plotspeclist, loc="upper left")


############################################################
# Prepare the species thermal conductivity data for plotting
# ==========================================================

# create arrays
# specify conductivity
kappa = np.zeros(points, dtype=np.double)
# start plotting loop #2
k = 0
# loop over selected gas species
for s in plotspeclist:
    # starting temperature at 300 [K]
    temp = 300.0
    # loop over temperature data points
    for i in range(points):
        conductivity = MyGasMech.species_cond(temp)
        index = MyGasMech.get_specindex(s)
        t[i] = temp
        # convert ergs to joules
        kappa[i] = conductivity[index] / ck.ERGS_PER_JOULE
        temp += dtemp
    plt.subplot(122)
    plt.plot(t, kappa, curvelist[k])
    k += 1
# plot conductivity versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Conductivity [J/cm-K-sec]")
plt.legend(plotspeclist, loc="upper left")


##########################################################################
# Evaluate the binary diffusion coefficients between different gas species
# ========================================================================
# Use the ``SpeciesDiffusionCoeffs()`` method to calculate
# the binary diffusion coefficients between pairs of gas species. Here,
# the binary diffusion coefficients are evaluated at 2 [atm] and 500 [K].
# The binary diffusion coefficient between CH\ :sub:`4` and O\ :sub:`2`
# is shown.

diffcoef = MyGasMech.species_diffusioncoeffs(2.0 * ck.P_ATM, 500.0)
id1 = MyGasMech.get_specindex(plotspeclist[0])
id2 = MyGasMech.get_specindex(plotspeclist[1])
c = diffcoef[id1][id2]
print(
    f"Diffusion coefficient for {plotspeclist[0]}"
    "against {plotspeclist[1]} is {c:e} [cm²/sec]."
)

# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_species_properties.png", bbox_inches="tight")
