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

"""Test for the mixture mixing calculations."""

from pathlib import Path

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# transport data not needed
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
# create the fuel mixture
# note: mixture pressures are not specified because pressure is not required
# for the calculations here
# the mixing process is assumed to take place at fixed pressure;
# i.e., the mixtures are at the same pressure
fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.x = [("CH4", 1.0)]
fuel.temperature = 300.0
# create the air mixture
air = ck.Mixture(MyGasMech)
# set mole fraction
air.x = [("O2", 0.21), ("N2", 0.79)]
air.temperature = 300.0
# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?)
mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermal_mixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)
# list molar composition
premixed.list_composition(mode="mole")
print()
# now create an argon mixture
ar = ck.Mixture(MyGasMech)
# species composition
ar.x = [("AR", 1.0)]
# mixture temperature
ar.temperature = 600.0
# dilute the premixed mixture adiabatically with the ar mixture by 30% by volume
dilute_recipe = [(premixed, 0.7), (ar, 0.3)]
# create the diluted mixture
diluted = ck.adiabatic_mixing(recipe=dilute_recipe, mode="mole")
# list molar composition
diluted.list_composition(mode="mole")
# show the mixture temperatures
print(f"the diluted mixture temperature is  {diluted.temperature:f} [K]")
print(f"the ar mixture temperature is       {ar.temperature:f} [K]")
print(f"the premixed mixture temperature is {premixed.temperature:f} [K]")

# return results for comparisons
resultfile = Path(current_dir) / "mixturemixing.result"
results = {}
results["state-temperature"] = [
    premixed.temperature,
    ar.temperature,
    float(diluted.temperature),
]
results["species-premixed_mole_fraction"] = premixed.x.tolist()
results["species-diluted_mole_fraction"] = diluted.x.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()
