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
.. _ref_create_mixture:

================
Create a mixture
================

A *mixture* is a core component of the PyChemkin framework. In addition to getting mixture
thermodynamic and transport properties, such as density, heat capacity, and viscosity, you can combine
two mixtures, find the equilibrium state of a mixture, or use a mixture to define the initial state of
a reactor. A PyChemkin *reactor model* is a black box that transforms a mixture
from its initial state to a new one.

The following schematic shows the basic operations available for a mixture in PyChemkin: **create**, **combine/mix**, and **transform** (by a reactor model).

.. figure:: mixture_concept.png
  :scale: 45 %
  :alt: the mixture concept

This example shows different ways to create a mixture in PyChemkin. The use of the composition *recipe* lets you provide just the non-zero species components with a list of species-fraction pairing tuples. Alternatively, the NumPy array lets you use a *full-size* (equal to the number of species) mole/mass fraction array to specify the mixture composition. The
``equivalence ratio()`` method creates a new mixture from predefined fuel and oxidizer mixtures by assigning an equivalence ratio value.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_create_mixture.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import copy
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
# The mechanism to load is the C2 NOx mechanism. This mechanism file and the associated data files
# come with the standard Ansys Chemkin installation in the ``/reaction/data`` directory.
# The C2 NOx mechanism file, in addition to the reactions, contains the thermodynamic
# and transport data of all species in the mechanism. In this case, you only need to specify
# the mechanism file, that is, ``chemfile``. If your simulation requires the transport properties, you
# must use the ``preprocess_transportdata()`` method to tell the PyChemkin preprocessor to also include the
# transport data. The preprocessor does not include the transport data by default.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the C2 NOx mechanism
MyGasMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files
# this mechanism file contains all the necessary thermodynamic and transport data
# thus, there is no need to specify thermodynamic and the transport data files
MyGasMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

# direct the preprocessor to include the transport properties (when the tran data file is not provided)
MyGasMech.preprocess_transportdata()

#####################################
# Preprocess the C2 NOx chemistry set
# ===================================
# The C2 NOx mechanism includes information about the *Soave* cubic
# Equation of State (EOS) for real-gas applications. The preprocessor
# indicates the availability of the real-gas model in the chemistry set processed.
# For example, during preprocessing, you see this print out: ``real-gas cubic EOS 'Soave' is
# available``.

# preprocess the mechanism files
iError = MyGasMech.preprocess()

######################################################################
# Set up gas mixtures based on the species in the C2 NOx chemistry set
# ====================================================================
# Create a gas mixture instance named ``premixed`` based on the ``My2ndMech`` chemistry set.

premixed = ck.Mixture(MyGasMech)
# set mixture pressure [dynes/cm2]
mixpressure = 2.0  # given in atm
# convert to dynes/cm2
premixed.pressure = mixpressure * ck.Patm
# set mixture temperature [K]
premixed.temperature = 500.0
# create a recipe for the molar composition of the mixture
mixture_recipe = [("CH4", 0.08), ("N2", 0.6), ("O2", 0.2), ("H2O", 0.12)]
# set mixture mole fractions
premixed.X = mixture_recipe

##################################
# Find mixture mean molecular mass
# ================================
# Use the ``WTM()`` method to get the mean molar mass of the gas mixture.

print(f"Mean molecular mass = {premixed.WTM:f} gm/mole")
print("=" * 40.0)

##############################
# List the mixture composition
# ============================
# Use the ``Y()`` method to automatically convert the mole fractions to mass fractions
# and vice versa. Use the ``list_composition()`` method to display
# only the non-zero components of the gas mixture.

print("mixture mass fractions (raw data):")
print(str(premixed.Y))
# switch back to mole fractions
print("\nmixture mole fractions (raw data):")
print(str(premixed.X))
# beautify the composition list
print("\nformatted mixture composition output:")
print("=" * 40)
premixed.list_composition(mode="mole")
print("=" * 40)

#################################
# Create a hard copy of a mixture
# ===============================
# Create a hard copy of the premixed mixture. (Soft copying, that is,
# using ``anotherpremixed = premixed``, also works.)
anotherpremixed = copy.deepcopy(premixed)
# display the molar composition of the new mixture
print("\nFormatted mixture composition of the copied mixture")
anotherpremixed.list_composition(mode="mole")
print("=" * 40.0)

#####################################
# Compute and plot mixture properties
# ===================================
# The PyChemkin ``Mixture`` module offers many basic methods to compute the
# thermodynamic and transport properties of a species and a gas mixture.
# The following code plots selected mixture properties as a function of the mixture
# temperature. It uses the ``RHO()`` and ``HML()`` methods to get the mixture density
# and mixture enthalpy, respectively. Use the ``mixture_viscosity()`` method to get
# the mixture transport property, ``viscosity``. Use the ``mixture_diffusion_coeffs()``
# method to get the mixture-averaged diffusion coefficient of CH\ :sub:`4`.
# The temperature and pressure are required to compute the properties.

plt.subplots(2, 2, sharex="col", figsize=(9, 9))
dTemp = 20.0
points = 100
# curve attributes
curvelist = ["g", "b--", "r:"]
# list of pressure values for the plot
press = [1.0, 5.0, 10.0]  # given in atm
# create arrays for the plot
# mixture density
rho = np.zeros(points, dtype=np.double)
# mixture enthalpy
enthalpy = np.zeros_like(rho, dtype=np.double)
# mixture viscosity
visc = np.zeros_like(rho, dtype=np.double)
# mixture averaged diffusion coefficient of CH4
diff_CH4 = np.zeros_like(rho, dtype=np.double)
# species index of CH4
CH4_index = MyGasMech.get_specindex("CH4")
# temperature data
T = np.zeros_like(rho, dtype=np.double)
# start the plotting of loop #1
k = 0
# loop over the pressure values
for j in range(len(press)):
    # starting temperature at 300K
    temp = 300.0
    # loop over temperature data points
    for i in range(points):
        # set mixture pressure [dynes/cm2]
        premixed.pressure = press[j] * ck.Patm
        # set mixture temperature [K]
        premixed.temperature = temp
        # get mixture density [gm/cm3]
        rho[i] = premixed.RHO
        # get mixture enthalpy [ergs/mol] and convert it to [kJ/mol]
        enthalpy[i] = premixed.HML() * 1.0e-3 / ck.ergs_per_joule
        # get mixture viscosity [gm/cm-sec]
        visc[i] = premixed.mixture_viscosity()
        # get mixture-averaged diffusion coefficient of CH4 [cm2/sec]
        diffcoeffs = premixed.mixture_diffusion_coeffs()
        diff_CH4[i] = diffcoeffs[CH4_index]
        T[i] = temp
        temp += dTemp

    # create subplots
    # plot mixture density versus temperature
    plt.subplot(221)
    plt.plot(T, rho, curvelist[k])
    plt.ylabel("Density [g/cm3]")
    plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper right")
    # plot mixture enthalpy versus temperature
    plt.subplot(222)
    plt.plot(T, enthalpy, curvelist[k])
    plt.ylabel("Enthalpy [kJ/mole]")
    plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
    # plot mixture viscosity versus temperature
    plt.subplot(223)
    plt.plot(T, visc, curvelist[k])
    plt.xlabel("Temperature [K]")
    plt.ylabel("Viscosity [g/cm-sec]")
    plt.legend(("1 atm", "5 atm", "10 atm"), loc="lower right")
    # plot mixture averaged CH4 diffusion coefficient versus temperature
    plt.subplot(224)
    plt.plot(T, diff_CH4, curvelist[k])
    plt.xlabel("Temperature [K]")
    plt.ylabel(r"$D_{CH_4}$ [cm2/sec]")
    plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
    k += 1

# plot legends
plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_create_mixture.png", bbox_inches="tight")
