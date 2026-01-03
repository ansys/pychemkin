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

r""".. _ref_heating_values:

=============================================
Calculate the heating values of fuel mixtures
=============================================

One of the advantages of PyChemkin is flexibility. You can establish your own
workflow with PyChemkin to meet your simulation goals. This tutorial is an example
of building a specialized algorithm to calculate the *heating values*,
the lower heating value (LHV) and the higher heating value (HHV), of fuel mixtures.

The **heating value** is the net heat release from the complete combustion of
a hydrocarbon (HC) in oxygen at the standard state condition (298.15 [K] and 1 [atm]).
The combustion products are mainly CO\ :sub:`2` and H\ :sub:`2`\ O, and are assumed
to be cooled back to the standard state condition in the heating value calculation.
The difference between the lower heating value and the higher heating value is
the final form of the H\ :sub:`2`\ O; the water is in *gas phase* for
the *lower* heating value, and in *liquid phase* for the *higher* heating value.
You can see that the higher heating value can be obtained by adding
the water *heat of vaporization* to the lower heating value. Thus,
the workflow for calculating the heating values of any HC fuels consists of
two main steps:

    1. Calculating the water heat of vaporization at the standard state condition:
    this can be done by getting the value from a well trusted database or
    by using the **chemkin** trick described below.

    2. Calculating the heat of combustion of the fuel mixture in pure oxygen:
    here you will use the ``find_equilibrium`` method with the *fixed pressure* and
    *fixed temperature* option (the default setting).

The lower heating value is the enthalpy different between the fresh fuel and
oxygen mixture and the final product mixture. Adding the heat of vaporization to
the lower heating value, and you get the higher heating value.

In this tutorial, you will compute the heating values of some pure fuel species
such as methane and n-butane as well as some fuel mixtures such as PRF RON 80 and
biodiesel. You can compare the values you get here with the known values from
a trusty database.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_heating_values.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

from pathlib import Path

import numpy as np  # number crunching

import ansys.chemkin.core as ck
from ansys.chemkin.core import Color
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.utilities import find_file

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
logger.debug("data directory: " + str(data_dir))
# set pressure & temperature condition of the standard state
thispressure = ck.P_ATM
thistemperature = 298.15

####################################
# Calculate the heat of vaporization
# ==================================
# Because the thermodynamic data (mainly the enthalpy) of both the vapor and the
# liquid water are available in the standard thermodynamic data file ``therm.dat``
# that comes with the Ansys Chemkin in the *reaction/data* directory, you can
#  in theory, compute the water heat of vaporization at a given temperature
# by finding the enthalpy difference between the water vapor and
# its liquid counterpart. In the ``getwaterheatofvaporization`` method,
# a mechanism with just the water vapor ``H2O`` and the liquid water ``H2O(L)``
# is created in situ. Simply find and return the enthalpy difference of
# the two species by using the ``SpecisH`` method after preprocessing the
# ``WaterMech`` ``Chemistry Set``.


def getwaterheatofvaporization(temp: float) -> float:
    """Compute water heat of vaporization [erg/g-water] at the given temperature."""
    """
    Use the enthalpy difference between water vapor and liquid water
    at the temperature. Enthalpy data depend on temperature only
    There are empirical formulas for heat of vaporization, for example,
    DIPPR EQ.

    Parameters
    ----------
        temp: double scalar
            water temperature [K]

    Returns
    -------
        enthalpy: double
            water enthalpy of vaporization [erg/g-water]

    """
    # compute water heat of vaporization
    # create a chemistry set object
    watermech = ck.Chemistry(label="Water Only")
    #
    # create a new mechanism input file
    #
    global current_dir
    waterfile = Path(current_dir) / "water_chem.inp"
    w = Path.open(waterfile, "w")
    # the water mechanism contains two species:
    # water vapor (H2O) and liquid water (H2O(L))
    # declare elements
    w.write("ELEMENT H O END\n")
    w.write("SPECIES\n")
    w.write("H2O   H2O(L)\n")
    w.write("END\n")
    # no reaction is needed for thermodynamic property calculation
    w.write("REACTION\n")
    w.write("END\n")
    # close the mechnaism file
    w.close()
    # set mechanism input files
    # including the full file path is recommended
    watermech.chemfile = str(waterfile)
    watermech.thermfile = str(data_dir / "therm.dat")
    # pre-process
    ierror = watermech.preprocess()
    if ierror != 0:
        return 0.0
    # get species enthalpies [erg/mol] at 298.15 [K]
    waterenthalpies = watermech.species_h(thistemperature)
    # compute heat of vaporization of water [erg/g-water]
    # H_water_vapor - H_liquid_water
    heatvaporization = (waterenthalpies[0] - waterenthalpies[1]) / watermech.wt[0]
    # remove the temporary water mechanism file
    Path(waterfile).unlink()
    return heatvaporization


##############################
# Create your own fuel library
# ==============================
# You can create a fuel "library" mechanism that contains typical fuel species,
# oxygen, and the complete combustion products o(carbon dioxide and water *vapor*)
# only. No reaction is needed for this purpose as you will perform
# equilibrium calculation to get the complete combustion product mixture.
#
# .. warning::
#   It is recommended **not** to include any intermediate species such as CH3
#   or OH in the fuel library as these species might interfere with
#   the equilibrium calculation used to determine the complete combustion products.
#
# .. note ::
#   The thermodynamic data file ``Gasoline-Diesel-Biodiesel_PAH_NOx_therm_MFL2024.dat``
#   should have the property data of most commonly seen fuel species.
#   This encrypted data file is developed under the *Model Fuel Library* project and
#   is available in the *reaction/data/ModelFuelLibrary*.
#

############################################
# Instantiate the **fuel** ``Chemistry Set``
# ==========================================
# Create the chemistry set object ``MyGasMech``. The mechanism input file
# ``fuel_chem.inp`` will be created below.
#
# .. note ::
#   You can keep this file for later uses, for example,
#   by adding more fuel species. (remember to remove the``remove(mymechfile))``
#   command at the end of this project).
#
MyGasMech = ck.Chemistry(label="EQ")

####################################
# Create the **fuel** mechanism file
# ==================================
#
# create a new mechanism input file
#
mymechfile = Path(current_dir) / "fuels_chem.inp"
m = Path.open(mymechfile, "w")
# the mechanism contains only the necessary species
# (fuel, oxygen, and major combustion products)
# declare elements
m.write("ELEMENT c h o END\n")
# declare species
# ch4: Methane                  c4h10: n-Butane
# nc5h12: n-Pentane             nc7h16: n-Heptane
# ic8h18: iso-Octane            nc9h20: n-Nonane
# nc10h22: n-Decane             hmn: Heptamethylnonane
# c6h5ch3: Toluene              c6h5c2h5: Ethylbenzene
# chx: Cyclohexane              mch: Methylcyclohexane
# decalin: Decalin              ch3oh: Methanol
# c2h5oh: Alcohol               ch3och3: Dimethyl ether (DME)
# nc4h9oh: n-Butanol            mb: Mythyl butanoate
# md: Methyl decanoate          mhd: Methyl stearate
m.write("SPECIES\n")
m.write("ch4 c4h10 nc5h12 nc7h16 ic8h18 nc9h20 nc10h22\n")
m.write("hmn c6h5ch3 c6h5c2h5 chx mch decalin\n")
m.write("ch3oh c2h5oh ch3och3 nc4h9oh\n")
m.write("mb md mhd\n")
m.write("o2 co2 h2o\n")
m.write("END\n")
# no reaction is needed for equilibrium calculation
m.write("REACTION\n")
m.write("END\n")
# close the mechnaism file
m.close()
#
# set mechanism input files
# including the full file path is recommended
# note that the "year" in thermodynamic data file name could be different.
# it's MFL2023 for Ansys Chemkin 2025R1, MFL2024 for 2025R2, ...
MyGasMech.chemfile = str(mymechfile)
therm_dir = str(data_dir / "ModelFuelLibrary" / "Full")
MyGasMech.thermfile = find_file(
    therm_dir,
    "Gasoline-Diesel-Biodiesel_PAH_NOx_therm_MFL",
    "dat",
)

############################################
# Pre-process the **fuel** ``Chemistry Set``
# ==========================================
# preprocess the fuel mechanism "MyGasMech`` just created.
ierror = MyGasMech.preprocess()
if ierror == 0:
    print(Color.GREEN + ">>> preprocess OK", end=Color.END)
else:
    print(Color.RED + ">>> preprocess failed!", end=Color.END)
    exit()

#####################################
# Find the water heat of vaporization
# ===================================
# Set the fresh fuel-oygen mixture pressure & temperature to
# the standard state condition for the heating value calculations.
# Since the difference between the *lower heating value (LHV)* and the
# *higher heating value (HHV)* is the final form of the water.
# You will calculate the LHV first and get the HHV by adding
# the water heat of vaporization to the LHV. Use the ``get_specindex``
# method to find the species index of water ``MyGasMech``.
#
# .. note ::
#   The water heat of vaporization is computed by the ``getwaterheatofvaporization``
#   method you created earlier in this project. Reactivate ``MyGasMech``
#   (``MyGasMech.activate``) after calling
#   ``getwaterheatofvaporization`` because another ``Chemistry Set`` ``WaterMech``
#   is used there.
#

# find the index for water vapor
watervapor_index = MyGasMech.get_specindex("h2o")

# water heat of vaporization [erg/g-water] at 298.15 [K]
# either call this method before creating the current Chemistry Set
# or use the activate method to switch back to the current Chemistry Set
# after the call
heatvaporization = getwaterheatofvaporization(thistemperature)

# switch back to MyGasMech
MyGasMech.activate()

#############################################
# Set up the unburned **fuel-oxygen** mixture
# ===========================================
# Set up the *unburned* "fuel-oxygen" mixture for
# the heating value calculations. Here the ``x_by_equivalence_ratio`` method
# with \ :math:`\phi = 1` is employed to generate a stoichiometric
# fuel-oxygen mixture. The advantage of using this method is that you do not
# need to calculate the "fuel-oxygen" composition for every fuel mixture.
# You can use a *lean* "fuel-oxygen" mixture, too, because the standard
# heat of formation of O\ :sub:`2` is zero by definition.
#
# Here you will compute the heating values of four fuel mixtures: CH\ :sub:`4`\ ,
# C\ :sub:`4`\ H\ :sub:`10`\ ,PRF 80
# (20% C\ :sub:`7`\ H\ :sub:`16` + 80% C\ :sub:`8`\ H\ :sub:`18`\ ), and a mock-up
# biodiesel mixture
# (90% C\ :sub:`19`\ H\ :sub:`38`\ O\ :sub:`2` + 10% CH\ :sub:`3`\ OH).

# prepare the fuel mixtures
fuel = ck.Mixture(MyGasMech)
fuel.pressure = thispressure
fuel.temperature = thistemperature
# list of fuel compositions (mole/volume fractions) of which the heating values
# will be computed.
# [Methane, n-Butane, PRF RON 80, biodiesel]
fuels = [
    [("ch4", 1.0)],
    [("c4h10", 1.0)],
    [("nc7h16", 0.2), ("ic8h18", 0.8)],
    [("mhd", 0.9), ("ch3oh", 0.1)],
]

# specify oxidizers = pure oxygen
oxid = ck.Mixture(MyGasMech)
oxid.x = [("o2", 1.0)]
oxid.pressure = thispressure
oxid.temperature = thistemperature
# get o2 index
oxy_index = MyGasMech.get_specindex("o2")

# specify the complete combustion product species
products = ["co2", "h2o"]
# no added species
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)

# instantiate the unburned fuel-oxygen mixture
unburned = ck.Mixture(MyGasMech)
unburned.pressure = thispressure
unburned.temperature = thistemperature

#####################################################
# Compute the fuel heating value of the fuel mixtures
# ===================================================
# Now compute the heating values: LHV and HHV, from the enthalpy different
# between the unburned stoichiometric fuel-oxygen mixture ``unburned``
# and the complete combustion product mixture ``burned``. The complete combustion
# product mixture is found by using the ``find_equilibrium`` method with
# the default *fixed pressure* and *fixed temperature* option, that is, the product
# mixture is also at the standard state condition.
lhv = np.zeros(len(fuels), dtype=np.double)
hhv = np.zeros_like(lhv, dtype=np.double)
fuelcount = 0
for f in fuels:
    # re-set the fuel composition (mole/volume fractions)
    fuel.x = f
    # create a soichiometric fuel-oxygen mixture
    ierror = unburned.x_by_equivalence_ratio(
        MyGasMech, fuel.x, oxid.x, add_frac, products, equivalenceratio=1.0
    )
    # get the mixture enthalpy of the initial mixture [erg/g]
    h_unburned = unburned.hml() / unburned.wtm

    # compute the complete combustion state (fixed temperature and pressure)
    # this step mimics the complete burning of the initial fuel-oxygen mixture
    # at constant pressure and the subsequent cooling of the combustion products
    # back to the original temperature
    burned = unburned.find_equilibrium()
    # get the mixture enthalpy of the final mixture [erg/g]
    h_burned = burned.hml() / burned.wtm

    # get total fuel mass fraction
    fmass = 0.0e0
    bmassfrac = unburned.y
    for i in range(MyGasMech.kk):
        if i != oxy_index:
            fmass += bmassfrac[i]
    # water vapor mass fraction in the urned mixture
    wmass = burned.y[watervapor_index]
    #
    if np.isclose(fmass, 0.0, atol=1.0e-10):
        # no fuel species exists in the unburned mixture
        print(f">>> error no fuel species in the unburned mixture {f}")
        exit()

    # compute the heating values [erg/g-fuel]
    lhv[fuelcount] = -(h_burned - h_unburned) / fmass
    hhv[fuelcount] = -(h_burned - (h_unburned + heatvaporization * wmass)) / fmass
    fuelcount += 1

############################
# Display the heating values
# ==========================
# List the fuel mixtures and their LHV and HHV. The heating values are converted
# from the cgs units [erg/g] to [kJ/g].
print(
    f"Fuel Heating Values at {thistemperature} [K] and {thispressure * 1.0e-6} [bar]\n"
)
for i in range(len(fuels)):
    print(f"fuel composition:  {fuels[i]}")
    print(f" LHV [kJ/g-fuel]:  {lhv[i] / ck.ERGS_PER_JOULE / 1.0e3}")
    print(f" HHV [kJ/g-fuel]:  {hhv[i] / ck.ERGS_PER_JOULE / 1.0e3}\n")

##########
# Clean up
# ========
# Delete arrays, mixture objects, and temporary files no longer needed.
del hhv, lhv
del fuel, oxid, unburned, burned
# delete the local mechanism file just created
Path(mymechfile).unlink()
