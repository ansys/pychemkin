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

"""Test for gas-phase reaction rate calculation."""

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
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
MyGasMech.tranfile = str(mechanism_dir / "grimech30_transport.dat")
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.x = [("CH4", 1.0)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 5.0 * ck.P_ATM
fuelmixture.temperature = 1500.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 5.0 * ck.P_ATM
air.temperature = 1500.0
# create the premixed mixture to be defined
premixed = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
ierror = premixed.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=1.0
)
if ierror != 0:
    raise RuntimeError
# list the composition of the premixed mixture
premixed.list_composition(mode="mole")
#
# compute reaction rates
#
# temperature and pressure are requires to compute the reaction rates
premixed.pressure = 5.0 * ck.P_ATM
premixed.temperature = 1600.0
# get the net species molar rates of production [mole/cm3-sec]
rop = premixed.rop()
# list the nonzero rates in descending order
print()
specrate_order, species_rates = premixed.list_rop()
# get the forward and the reverse rates of each reaction
kf, kr = premixed.rxn_rates()
print()
print(f"reverse reaction rates: (raw values of all {MyGasMech.ii_gas:d} reactions)")
print(str(kr))
print("=" * 40)
# list the nonzero net reaction rates
rxn_order, net_rxn_rates = premixed.list_reaction_rates()
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
# change the mixture temperature
premixed.temperature = 1800.0
# get the list the nonzero net reaction rates at the new temperature
rxn_order, net_rxn_rates = premixed.list_reaction_rates()
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
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("reaction_rates.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "reactionrates.result"
results = {}
results["state-order_1800"] = rxn_order.tolist()
results["rate-net_reaction_rate_1800"] = net_rxn_rates.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()
