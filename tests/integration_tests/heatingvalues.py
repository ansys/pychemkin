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

"""Fuel species heating value test for the equilibrium calculation."""

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
global data_dir
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
logger.debug("data directory: " + str(data_dir))
# set pressure & temperature condition
thispressure = ck.P_ATM
thistemperature = 298.15


#
def getwaterheatofvaporization(temp) -> float:
    """Compute water heat of vaporization [erg/g-water] at the given temperature."""
    """
    Use the enthalpy difference between water vapor and liquid water at the temperature
    Enthalpy data depend on temperature only
    There are empirical formulas for heat of vaporization, for example, DIPPR EQ

    Parameters
    ----------
        temp: double scalar
            water temperature [K]

    Returns
    -------
        heatvaporization: double
            heat of vaporrization of water at the given temperature [ergs/g]
    """
    # compute water heat of vaporization
    # create a chemistry set object
    watermech = ck.Chemistry(label="Water Only")
    #
    # create a new mechanism input file
    #
    waterfile = Path(current_dir) / "water_chem.inp"
    w = waterfile.open(mode="w")
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
    waterenthalpies = watermech.species_h(temp)
    # compute heat of vaporization of water [erg/g-water]
    # H_water_vapor - H_liquid_water
    heatvaporization = (waterenthalpies[0] - waterenthalpies[1]) / watermech.wt[0]
    # remove the temporary water mechanism file
    Path(waterfile).unlink()
    return heatvaporization


#
# create the mechanism file with fuel species and complete combustion products only
# no reaction
# create a chemistry set object
MyGasMech = ck.Chemistry(label="EQ")
#
# create a new mechanism input file
#
mymechfile = Path(current_dir) / "fuels_chem.inp"
m = mymechfile.open(mode="w")
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
MyGasMech.chemfile = str(mymechfile)
therm_dir = str(data_dir / "ModelFuelLibrary" / "Full")
MyGasMech.thermfile = find_file(
    therm_dir,
    "Gasoline-Diesel-Biodiesel_PAH_NOx_therm_MFL",
    "dat",
)

# pre-process
ierror = MyGasMech.preprocess()
if ierror == 0:
    print(Color.GREEN + ">>> preprocess OK", end=Color.END)
else:
    print(Color.RED + ">>> preprocess failed!", end=Color.END)
    exit()
#
exit()
# set pressure & temperature condition
thispressure = ck.P_ATM
thistemperature = 298.15
# create the unburned fuel-oxygen mixture
unburned = ck.Mixture(MyGasMech)
unburned.pressure = thispressure
unburned.temperature = thistemperature
# find the index for water vapor
watervapor_index = MyGasMech.get_specindex("h2o")
# water heat of vaporization [erg/g-water] at 298.15 [K]
# either call this method before creating the current Chemistry Set
# or use the activate method to switch back to the current Chemistry Set after the call
heatvaporization = getwaterheatofvaporization(thistemperature)
# switch back to MyGasMech
MyGasMech.activate()
# prepare the fuel mixture
fuel = ck.Mixture(MyGasMech)
fuel.pressure = thispressure
fuel.temperature = thistemperature
# list of fuel compositions (mole/volume fractions)
# of which the heating values will be computed
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
#
# compute fuel heating values
#
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
    bmassfrac = unburned.Y
    for i in range(MyGasMech.kk):
        if i != oxy_index:
            fmass += bmassfrac[i]
    # water vapor mass fraction in the urned mixture
    wmass = burned.Y[watervapor_index]
    #
    if np.isclose(fmass, 0.0, atol=1.0e-10):
        # no fuel species exists in the unburned mixture
        print(f">>> error no fuel species in the unburned mixture {f}")
        exit()
    # compute the heating values [erg/g-fuel]
    lhv[fuelcount] = -(h_burned - h_unburned) / fmass
    hhv[fuelcount] = -(h_burned - (h_unburned + heatvaporization * wmass)) / fmass
    fuelcount += 1

# display results
print(
    f"Fuel Heating Values at {thistemperature} [K] and {thispressure * 1.0e-6} [bar]\n"
)
for i in range(len(fuels)):
    print(f"fuel composition:  {fuels[i]}")
    print(f" LHV [kJ/g-fuel]:  {lhv[i] / ck.ERGS_PER_JOULE / 1.0e3}")
    print(f" HHV [kJ/g-fuel]:  {hhv[i] / ck.ERGS_PER_JOULE / 1.0e3}\n")

# return results for comparisons
resultfile = Path(current_dir) / "heatingvalues.result"
results = {}
lhv = lhv / ck.ERGS_PER_JOULE / 1.0e3
results["state-LHV"] = lhv.tolist()
hhv = hhv / ck.ERGS_PER_JOULE / 1.0e3
results["state-HHV"] = hhv.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()
