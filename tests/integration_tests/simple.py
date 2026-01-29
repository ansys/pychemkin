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

"""Test for defining and preprocessing gas phase reaction mechanism."""

from pathlib import Path

import ansys.chemkin.core  # import PyChemkin
from ansys.chemkin.core.logger import logger

# create a Chemistry Set for GRI 3.0 mechanism in the data directory
mechanism_dir = Path(ansys.chemkin.core.ansys_dir) / "reaction" / "data"
# set up mechanism file names
mech_file = str(mechanism_dir / "grimech30_chem.inp")
therm_file = str(mechanism_dir / "grimech30_thermo.dat")
tran_file = str(mechanism_dir / "grimech30_transport.dat")
# instantiate Chenistry Set 'GasMech'
GasMech = ansys.chemkin.core.Chemistry(
    chem=mech_file, therm=therm_file, tran=tran_file, label="GRI 3.0"
)
# pre-process the Chemistry Set
status = GasMech.preprocess()
# check preprocess status
if status != 0:
    # failed
    print(f"Preprocessing: Error encountered. Code = {status:d}.")
    print(f"see the summary file {GasMech.summaryfile} for details")
    logger.error("PreProcess failed")
    exit()
# Create Mixture 'air' based on 'GasMech'
air = ansys.chemkin.core.Mixture(GasMech)
# set 'air' condition
# mixture pressure in [dynes/cm2]
air.pressure = 1.0 * ansys.chemkin.core.P_ATM
# mixture temperature in [K]
air.temperature = 300.0
# mixture composition in mole fractions
air.x = [("O2", 0.21), ("N2", 0.79)]
#
print(f"pressure    = {air.pressure / ansys.chemkin.core.P_ATM} [atm]")
print(f"temperature = {air.temperature} [K]")
# print the 'air' composition in mass fractions
air.list_composition(mode="mass")
# get 'air' mixture density [g/cm3]
print(f"the mixture density   = {air.rho} [g/cm3]")
# get 'air' mixture viscosity [g/cm-sec] or [poise]
print(f"the mixture viscosity = {air.mixture_viscosity() * 100.0} [cP]")

# return results for comparisons
current_dir = Path.cwd()
resultfile = current_dir / "simple.result"
results = {}
results["state-temperature"] = [air.temperature]
results["state-pressure"] = [air.pressure]
results["state-density"] = [air.rho]
results["state-viscosity"] = [air.mixture_viscosity() * 100.0]
results["species-mole_fraction"] = air.x.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()
