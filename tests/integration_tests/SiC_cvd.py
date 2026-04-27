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

"""Simulation the Chemical Vapor Deposition (CVD) process."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.core.stirreactors.PSR import PSRSetVolumeFixedTemperature as Psr

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

# set numpy printing option
np.set_printoptions(legacy="1.25")

# Create a chemistry set
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the location of this example py file
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir) / ".." / "data"
script_dir = str(script_dir_obj)
# use the relative path to locate the example mechanism folder
example_mech_data = script_dir_obj / ".." / "data"
# create a chemistry set based on the Si-N CVD mechanism
MyCVDMech = ck.Chemistry(label="Si-N_CVD")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "SIN_CVD"
MyCVDMech.chemfile = str(example_mech_data / local_data_dir / "chem_SIN_cvd.inp")
MyCVDMech.surffile = str(example_mech_data / local_data_dir / "surf_SIN_cvd.inp")
MyCVDMech.thermfile = str(data_dir / "therm.dat")

# preprocess the mechanism files
_ = MyCVDMech.preprocess()

# Set up the precursors stream
# precursor stream compositions
inlet_recipes = [
    [("SIF4", 0.14286), ("NH3", 0.85714)],
    [("SIF4", 0.16), ("NH3", 0.84)],
    [("SIF4", 0.18), ("NH3", 0.82)],
    [("SIF4", 0.20), ("NH3", 0.80)],
    [("SIF4", 0.22), ("NH3", 0.78)],
]
# create the gas precursors mixture
precursors = Stream(MyCVDMech)
# set fuel composition: hydrogen diluted by nitrogen
# use the first inlet composition to instantiate the inlet stream object
precursors.x = inlet_recipes[0]
# setting the pressure and temperature is not required in this case
precursors.pressure = 1.8 * ck.P_TORRS
precursors.temperature = 1713.15  # inlet temperature
precursors.sccm = 1.13e4  # cm3/sec @ standard condition

# Create the PSR to predict the gas composition of the outlet stream
cvd_reactor = Psr(precursors, label="CVD_PSR")

# Set up additional reactor model parameters
# set PSR volume (cm3): required for PSRSetVolumeFixedTemperature model
cvd_reactor.volume = 2000.0

# Connect the inlet to the reactor
cvd_reactor.set_inlet(precursors)

# Set up surface chemistry parameters
# set PSR total active internal surface area (cm2):
cvd_reactor.area = 950.0
# get the number of surface material defined in the surface mechanism
n_material = cvd_reactor.get_numb_material
# get all surface material names
cvd_material_names = cvd_reactor.get_material_names()
# set surface area fraction (by default, materials have the same area)
area_frac = 1.0 / n_material
# site species symbols: list[list[str]]
site_names = cvd_reactor.get_site_species_names()
# bulk species symbols: list[list[str]]
bulk_names = cvd_reactor.get_bulk_species_names()
# list the material names and the site species and the bulk species symbols
# per material
print(f"number of surface materials = {n_material}")
for m in range(n_material):
    print(f"material names: {cvd_material_names[m]}")
    print(f"   list of site species:\n     {site_names[m]}")
    print(f"   list of bulk species:\n     {bulk_names[m]}")

# set up surface coverage of all surface materials
site_recipe = [[("HN_SIF(S)", 0.5), ("HN_NH2(S)", 0.5)]]
bulk_recipe = [[("SI(D)", 1.0), ("N(D)", 1.0)]]
# set the surface condition of the surface material
for i, mname in enumerate(cvd_material_names):
    # set area fraction [-] (optional: by default, materials have the same area)
    cvd_reactor.set_material_area_fraction(mname, area_frac)
    # set material temperature [K] if it is different from the gas temperature
    # (optional: by default, it is the same as the mixture temperature)
    cvd_reactor.set_material_temperature(mname, precursors.temperature)
    # set reactor surface coverage of the material
    # site fractions (optional: by default, all sites have the same fraction)
    cvd_reactor.set_site_fraction(mname, site_recipe[i])
    # bulk activities (optional: by default, all bulk species activities = 1.0)
    cvd_reactor.set_bulk_activity(mname, bulk_recipe[i])
# scale surface reaction rates
cvd_reactor.surface_ratemultiplier = 1.0

# Set steady-state solver controls
# reset the tolerances in the steady-state solver
cvd_reactor.steady_state_tolerances = (1.0e-12, 1.0e-8)
cvd_reactor.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
cvd_reactor.set_species_floor(-1.0e-10)

# Run the inlet precursor composition parameter study
# solution arrays
numb_runs = len(inlet_recipes)
# get gas-species index
i_sif4 = precursors.get_specindex("SIF4")
i_nh3 = precursors.get_specindex("NH3")
i_hf = precursors.get_specindex("HF")
# get surface species index
i_hn_sif_global, i_hn_sif_local = precursors.get_surf_specindex("HN_SIF(S)")
i_f2sinh_global, i_f2sinh_local = precursors.get_surf_specindex("F2SINH(S)")
i_hn_fsinh_2_global, i_hn_fsinh_2_local = precursors.get_surf_specindex("HN(FSINH)2(S)")
i_si_d_global, i_si_d_local = precursors.get_surf_specindex("SI(D)")
# solution arrays
# SiF4 inlet mole fraction
sif4_inlet = np.zeros(numb_runs, dtype=np.double)
# staedy-state gas mole fraction of HF
hf_ss_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# steady-state site fraction of HN_SiF(S)
hn_sif_s_ss_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# steady-state site fraction of F2SiNH(S)
f2sinh_s_ss_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# steady-state site fraction of HN(FSiNH)2(S)
hn_fsinh_2_s_ss_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# net surface rate of production [mole/cm2-sec] of NH3
nh3_rop_solution = np.zeros_like(sif4_inlet, dtype=np.double)
sif4_rop_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# linear deposition rate [micro/min] of Si(D)
si_d_growth_solution = np.zeros_like(sif4_inlet, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet precursor compositions
for i, r in enumerate(inlet_recipes):
    # update inlet stream composition and the estimated reactor gas composition
    precursors.x = r
    cvd_reactor.reset_inlet(precursors)
    cvd_reactor.reset_estimate_composition(r)
    # reset estimated surface coverage
    for n, mname in enumerate(cvd_material_names):
        cvd_reactor.set_site_fraction(mname, site_recipe[n])
        cvd_reactor.set_bulk_activity(mname, bulk_recipe[n])
    # run the PSR model
    runstatus = cvd_reactor.run()
    # check run status
    if runstatus != 0:
        # Run failed.
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # Run succeeded.
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # postprocess the solution profiles
    solnmixture = cvd_reactor.process_solution()
    # print the steady-state solution values
    # print(f"Steady-state temperature = {solnmixture.temperature} [K].")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    sif4_inlet[i] = precursors.x[i_sif4]
    # gas species solution
    hf_ss_solution[i] = solnmixture.x[i_hf]
    # get surface properties (use local index)
    s_frac = solnmixture.surface_chemistry.get_all_site_frac()
    b_rate = solnmixture.surface_chemistry.get_all_bulk_growth_rates()
    # surface site species fraction
    hn_sif_s_ss_solution[i] = s_frac[i_hn_sif_local]
    f2sinh_s_ss_solution[i] = s_frac[i_f2sinh_local]
    hn_fsinh_2_s_ss_solution[i] = s_frac[i_hn_fsinh_2_local]
    # linear growth rate of bulk species [cm/sec] -> [micron/min]
    si_d_growth_solution[i] = b_rate[i_si_d_local] * 6.0e5
    # get gas species surface chemistry ROP [mole/cm2-sec]
    nh3_rop_solution[i] = 0.0
    sif4_rop_solution[i] = 0.0
    for m in cvd_material_names:
        surf_rop, _ = solnmixture.rop_surf(m)
        # gas species net rate of production due to surface chemistry
        nh3_rop_solution[i] += surf_rop[i_nh3]
        sif4_rop_solution[i] += surf_rop[i_sif4]

# compute the total runtime
runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec] over {numb_runs} runs.")

# clean up
ck.done()

# Plot the parameter study results
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("PSR CVD Reactor Solution", fontsize=16)
plt.subplot(221)
plt.plot(sif4_inlet, hf_ss_solution, "b-")
plt.ylabel("HF Mole Fraction")
plt.subplot(222)
plt.plot(sif4_inlet, hn_sif_s_ss_solution, "b-")
plt.plot(sif4_inlet, f2sinh_s_ss_solution, "r-")
plt.plot(sif4_inlet, hn_fsinh_2_s_ss_solution, "g-")
plt.legend(["HN_SiF(S)", "F2SiNH(S)", "HN(FSiNH)2(S)"], loc="lower right")
plt.yscale("log")
plt.ylabel("Site Fraction")
plt.subplot(223)
plt.plot(sif4_inlet, nh3_rop_solution, "g-")
plt.plot(sif4_inlet, sif4_rop_solution, "r-")
plt.legend(["NH3", "SIF4"], loc="lower left")
plt.xlabel("Inlet SiF4 Mole Fraction")
plt.ylabel("Surface ROP [mole/cm2-sec]")
plt.subplot(224)
plt.plot(sif4_inlet, si_d_growth_solution, "b-")
plt.xlabel("Inlet SiF4 Mole Fraction")
plt.ylabel("Si(D) Linear Growth Rate [micron/min]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_PSR_cvd.png", bbox_inches="tight")

resultfile = Path(current_dir) / "PSR_cvd.result"
results = {}
results["species-SiF4_inlet_mole_fraction"] = sif4_inlet.tolist()
results["species-HF_mole_fraction"] = hf_ss_solution.tolist()
results["species-HN_SiF(S)_site_fraction"] = hn_sif_s_ss_solution.tolist()
results["rate-NH3_rop"] = nh3_rop_solution.tolist()
results["rate-Si(D)_linear_growth_rate"] = si_d_growth_solution.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()
