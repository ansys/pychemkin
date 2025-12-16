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
.. _ref_reaction_rates:

===================
Rank reaction rates
===================
When you have a mixture in PyChemkin, you can not only obtain its properties but also extract the rate information. The mixture can be any one of the following:

- Set up from scratch.
- Obtained from certain mixture operations.
- Based on a point (time or grid) solution of a reactor simulation.

The mixture rate utilities let you access the net production rate of species (ROP), forward and reverse reaction rates per reaction (RR), and net chemical heat release rate (HRR).

This example shows how to use PyChemkin mixture rate tools to derive useful information
from the raw data, such as isolating dominant reactions at different mixture conditions.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_reaction_rates.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
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
# The first mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")


##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
iError = MyGasMech.preprocess()


################################################################
# Set up gas mixtures based on the species in this chemistry set
# ==============================================================
# Before you can create a premixed fuel-oxidizer mixture by the equivalence ratio, you must
# first create fuel and the oxidizer mixtures.


#########################
# Create the fuel mixture
# =======================
# The fuel mixture consists of 100% methane.

fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.X = [("CH4", 1.0)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 5.0 * ck.Patm
fuelmixture.temperature = 1500.0


########################
# Create the air mixture
# ======================
# The air mixture consists of oxygen and nitrogen. The temperature and the pressure of
# this mixture are set to the values of the fuel mixture.

air = ck.Mixture(MyGasMech)
air.X = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 5.0 * ck.Patm
air.temperature = 1500.0


##############################################
# Define the combustion products and additives
# ============================================
# To use the equivalence ratio method, you must define
# the complete combustion products and the composition
# of the additives.

# products from the complete combustion of the fuel and air mixtures
products = ["CO2", "H2O", "N2"]
# Specify mole fractions of the added/inert mixture. You can also create an additives mixture here.
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros


###############################################
# Create a placeholder for the fuel-air mixture
# =============================================
# You can create an empty mixture object that is fully defined later.
# The composition of the ``premixed`` mixture is determined when calling one of the
# mixture composition methods:
#
# - ``X()``
# - ``Y()``
# - ``X_byEquivalence_Ratio()``
# - ``Y_byEquivalence_Ratio()``
#
# When calling the ``X_byEquivalence_Ratio()`` or ``Y_byEquivalence_Ratio()`` method,
# the equivalence ratio is specified by this parameter: ``equivalenceratio=1.0``.

# create the premixed mixture to define
premixed = ck.Mixture(MyGasMech)
# define the actual composition by the equivalence ratio
iError = premixed.X_by_Equivalence_Ratio(
    MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=1.0
)
if iError != 0:
    # check fuel-air mixture creation status
    print("Error: Failed to create the fuel-air mixture.")
    exit()


#######################################################
# Display the molar composition of the premixed mixture
# =====================================================
# List the composition of the premixed mixture for verification.
premixed.list_composition(mode="mole")


#####################################################################
# Evaluate reaction rates at a given mixture temperature and pressure
# ===================================================================
# Compute and rank the *net* reaction rates in the ``MyGasMech`` chemistry set at 1600 [K]
# and 5 [atm]. The *net* rate of a reaction is computed by summing the forward rate ``kf``
# and the reverse rate ``kr``. The ``RxnRates`` rate utility returns both the forward and
# the reverse rates of each reaction in the mechanism.
#
# Get the net production rates of *all* species by using the ``ROP`` method.
# The ROP values obtained are in [mole/cm\ :sup:`3`\ -sec].
# Alternatively, you can use the ``massROP`` method to get mass-unit ROP values in
# [g/cm\ :sup:`3`\ -sec].
#
# Use the optional ``threshold`` parameter to get only the ROP values with absolute value above
# the given threshold value. For example, ``premixed.ROP(threshold=1.0e-8)`` returns a
# ``rop`` array in which any raw ROP value with an absolute value less than 1.0e-8 is set to zero.
# By default, ``threshold=0.0``, in which case all raw ROP values are returned.
#
# .. note::
#    Temperature and pressure are required to compute the reaction rates.

# set the temperature and the pressure of the premixed mixture
premixed.pressure = 5.0 * ck.Patm
premixed.temperature = 1600.0

########################################################
# Obtain and sort the rates of production of all species
# ======================================================
# Use the ``ROP()`` method to get the rates of production (ROP) of all species in the
# ``MyGasMech`` chemistry set. You can use the ``list_ROP()`` method to list and sort
# the ROP values. ``spec_rate_order`` contains the sorted species indices in descending
# order. ``species_rates`` contains the sorted ROP values.

rop = premixed.ROP()
# list the non-zero rates in descending order
print()
spec_rate_order, species_rates = premixed.list_ROP()


##########################################################
# Evaluate and sort the rates of production of all species
# ========================================================
# Use the ``RxnRates()`` method to get the forward rate (``kf``)
# and reverse rate (``kr``) of each reaction. The ``list_reaction_rates()``
# method computes the net rate of each reaction and sorts them in ascending
# or descending order. The net reaction rate is computed by summing the forward and
# reverse rates of each reaction. ``rxn_order`` is a list of the ranked reaction indices.
# ``net_rxn_rates`` contains the net reaction rates in the same order.
#
# .. note::
#   The species and the reactions indexes returned by the ``list_ROP()`` and ``list_reaction_rates()``
#   methods are 0-based indexes. Increment the returned index by 1 to get the actual 1-based index.

kf, kr = premixed.RxnRates()
print()
print(f"Reverse reaction rates: (raw values of all {MyGasMech.IIGas:d} reactions)")
print(str(kr))
print("=" * 40)
# list the non-zero net reaction rates
rxn_order, net_rxn_rates = premixed.list_reaction_rates()


############################################
# Plot the sorted reaction rates at 1600 [K]
# ==========================================
# Display the top net reaction rates and their reaction strings in the commonly used horizontal
# bar plot. Use the ``get_gas_reaction_string`` utility of the chemistry set to get the actual
# reaction for the Y-axis label.

# create a rate plot
plt.rcParams.update({"figure.autolayout": True})
plt.subplots(2, 1, sharex="col", figsize=(10, 5))
# convert reaction # from integers to strings
rxnstring = []
for i in range(len(rxn_order)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(rxn_order[i] + 1))
# use horizontal bar chart
plt.subplot(211)
plt.barh(rxnstring, net_rxn_rates, color="blue", height=0.4)
# use log scale on x axis
plt.xscale("symlog")
# plt.ylabel('reaction')
plt.text(-3.0e-4, 0.5, "T = 1600K", fontsize=10)


#####################################################################
# Evaluate reaction rates at a given mixture temperature and pressure
# ===================================================================
# Compute and rank the net reaction rates in the ``MyGasMech`` chemistry set at 1800 [K]
# and 5 [atm].

# change the mixture temperature
premixed.temperature = 1800.0


######################################################################
# Evaluate and sort the rates of production of all species at 1800 [K]
# ====================================================================
# Get the list of  non-zero net reaction rates at the new temperature.
rxn_order, net_rxn_rates = premixed.list_reaction_rates()


############################################
# Plot the sorted reaction rates at 1800 [K]
# ==========================================
# Display the top net reaction rates and their reaction strings in the commonly used horizontal bar plot.

plt.subplot(212)
# convert reaction # from integers to strings
rxnstring.clear()
for i in range(len(rxn_order)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(rxn_order[i] + 1))
plt.barh(rxnstring, net_rxn_rates, color="orange", height=0.4)
plt.xlabel("reaction rate [mole/cm3-sec]")
# plt.ylabel('reaction')
plt.text(-3.0e-4, 0.5, "T = 1800K", fontsize=10)
# use log scale on x axis
plt.xscale("symlog")

# plot both ranking results
if interactive:
    plt.show()
else:
    plt.savefig("plot_reaction_rates.png", bbox_inches="tight")
