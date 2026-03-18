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

r""".. _ref_PSR_CVD:

======================================================
Simulating the Chemical Vapor Deposition (CVD) process
======================================================
The chemical vapor deposition (CVD) process plays an important part in the
semiconductor industries, especially for manufacturing and processing the wafers.

This example project presents a model for the CVD of silicon nitride (Si-N) in a
steady-state PSR. The utilization of the simple PSR model is justified because
the CVD process is operated under the "kinetic limited" condition and is running
steadily. The CVD PSR is characterized by the volume, the substrate surface area,
and the inlet precursor gas volumetric flow rate; consequently the reactor
residence time remains constant even when the precursor stream composition or
the reactor temperature is changed.

The current CVD process operates at a low pressure of 1.8 Torr, and a high temperature
(for both gas phase and the substrate surface) of 1440 |deg|C. The precursor stream
consists of SiF\ :sub:`4` and NH\ :sub:`3`\ .
The chemistry set contains the gas phase mechanism input file, "chem_SIN_cvd.inp",
the surface mechanism input file "surf_SIN_cvd.inp", and the thermodynamic data file,
"therm.dat". The two mechanism input files are in the "example\data\SIN_CVD" folder,
and the generic thermodynamic data file is in the "data" folder of the Ansys Chemkin
installation. The chemistry set is described in the Si-N CVD from Silicon Tetrafluoride
and Ammonia section of the *Chemkin Tutorial* manual (p. 346).

When a reactor model is instantiated with a ``Mixture`` or a ``Stream`` that is made of
a ``Chemistry Set`` with a surface mechanism, the surface chemistry module of the
reactor model will be automatically initialized with the surface mechanism data from
the ``Chemistry Set``. There is no additional step to "turn on" the surface mechanism
for the reactor model.

You will determine the effects of the inlet SiF\ :sub:`4` and NH\ :sub:`3` molar
ratio on the linear growth rate (deposition rate) of the Si-N layers on
the substrate while keeping the reactor temperature and pressure fixed. The inlet
composition can be easily changed by resetting the inlet stream composition.

In addition to the bulk species (Si and N) growth rates, you can get the reactor
gas composition and the substrate surface coverage from the steady-state PSR solution
for each case. You can also obtain the total surface rate-of-production
(surface ROPs) from the solutions to analyze the important reaction pathways of the CVD
process.

The solution plots show the dependencies of the mole fraction of gas product HF, the
substrate surface coverage, the Si linear growth rate, and the surface consumption
rates of the precursors on the inlet SiF\ :sub:`4` molar ratios. The results shows that
the deposition rate can be controlled by adjusting the precursor molar ratio. The growth
rate could reach its max value When both precursors have the same molar concentration.
Moreover, the negative surface ROPs indicate that surface reactions (rather than gas
reactions) are the dominate decomposition pathways of the gas-phase precursors.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_PSR_cvd.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

from ansys.chemkin.core.inlet import Stream

from ansys.chemkin.core.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.core.stirreactors.PSR import PSRSetVolumeFixedTemperature as Psr
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

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
# This example uses the surface mechanism of Si-N chemical vapor deposition (Si-N CVD)
# from Silicon tetrafluoride (SiF\ :sub:`4`\ `). The reaction mechanism consists of two
# parts: the gas-phase mechanism, ``chem_CIN_cvd.inp``, describing the breakdown
# of the precursors in the gas phase; and the surface mechanism, ``surf_SIN_cvd.inp``,
# describing the reactions taking place on the substrate surface, such as the
# adsorption of the gas precursors, the reactions between the adsorbed surface species,
# and desorption of gas products for the overall process of the deposition of the Si
# and the N atoms.
#
# This mechanism is provided in the *examples\data* folder of the local
# *PyChemkin* installation. You must provide the mechanism file names,
# ``chem_SIN_cvd.inp`` and ``surf_SIN_cvd.inp``, with the correct path
# ``example_mech_data`` to the Chemistry Set before pre-processing. The thermodynamic
# data file from the default Ansys Chemkin installation will be used.
#
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the location of this example py file
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir) / ".." / "data"
script_dir = str(script_dir_obj)
# use the relative path to locate the example mechanism folder
example_mech_data = script_dir_obj / ".." / "data"
print(f"example data = {str(example_mech_data)}")
# create a chemistry set based on the Si-N CVD mechanism
MyCVDMech = ck.Chemistry(label="Si-N_CVD")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "SIN_CVD"
MyCVDMech.chemfile = str(example_mech_data / local_data_dir / "chem_SIN_cvd.inp")
MyCVDMech.surffile = str(example_mech_data / local_data_dir / "surf_SIN_cvd.inp")
MyCVDMech.thermfile = str(data_dir / "therm.dat")

#######################################
# Preprocess the CVD chemistry set
# =====================================
# preprocess the mechanism files
_ = MyCVDMech.preprocess()

#####################################
# Set up the precursors stream
# ===================================
# Instantiate a stream feed for the inlet gas mixture.
# The stream is a mixture with the addition of the
# inlet flow rate. You specify inlet gas properties in the same way as you
# set up a mixture.
#
# Use the ``mass_flowrate`` method to assign the inlet mass flow rate.
# In this example, the ``sccm`` method, [cm3/min] @ the standard condition,
# is used specify the inlet volumetric flow rate. And the different inlet
# precursors compositions are given by a list of mole fractions recipes
# ``inlet_recipes``.

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

####################################################################
# Create the PSR to predict the gas composition of the outlet stream
# ===================================================================
# Instantiate the PSR ``cvd_reactor`` as a ``PSRSetVolumeFixedTemperature``
# object. Since the goal is to show how precursor stream
# composition affects the Si-N deposition rate, both the gas and the surface
# temperatures are fixed. The inlet precursor mixture ``precursors``
# composition will be applied as the estimated reactor condition of the
# ``cvd_reactor`` by default. You can overwrite
# any estimated reactor conditions by using appropriate methods. For example,
# ``cvd_reactor.temperature = 2000.0`` changes the estimated reactor temperature
# from 1713.15 to 2000 [K]. The residence time is calculated from the inlet steam
# velocity and the PSR volume.

cvd_reactor = Psr(precursors, label="CVD_PSR")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For a steady-state PSR, you must provide
# either the residence time or the reactor volume. You can also make changes to
# any estimated reactor conditions, such as the reactor temperature, the gas
# composition, and the surface coverage, to improve the solver convergence performance
# in addition to the solver parameters.

# set PSR volume (cm3): required for ``PSRSetVolumeFixedTemperature`` model
cvd_reactor.volume = 2000.0

##################################
# Connect the inlet to the reactor
# ================================
# You must connect at least one inlet to the open reactor. Use the ``set_inlet()``
# method to add a stream to the PSR. Inversely, use the ``remove_inlet()`` to
# disconnect an inlet from the PSR.
#
# .. note ::
#   There is no limit on the number of inlets that can be connected to a PSR.
#

# connect the inlet to the reactor
cvd_reactor.set_inlet(precursors)

#####################################
# Set up surface chemistry parameters
# ===================================
# When the ``Mixture`` or the ``Stream`` used to initiate the reactor model includes
# *surface chemistry*, the surface chemistry data and properties will automatically
# set up by the reactor model. However, you need to provide some essential surface
# chemistry related model parameters. For example, you need to specify the activate
# surface area as well as the temperature for each material defined in the surface
# mechanism if it is different from the gas temperature. By default,
# all site species have the same fraction and all bulk have activity of unity.
# You can explicitly provide the site fractions and bulk activities when the default
# surface coverage fail to get a converged solution.

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

##################################
# Set steady-state solver controls
# ================================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances. Here the tolerances that the steady-state solver
# is to use for the *steady-state search* and the *pseudo time stepping* stages
# are changed by using the ``steady_state_tolerances`` method. Sometimes, during
# the iterations, some species mass fractions might become negative, causing the
# solver to report an error and stop. To overcome this issue, you can use the
# ``set_species_floor()`` method to provide a small cushion, allowing species
# mass fractions to go slightly negative (during the iterations) by resetting
# the mass fraction floor value.

# reset the tolerances in the steady-state solver
cvd_reactor.steady_state_tolerances = (1.0e-12, 1.0e-8)
cvd_reactor.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
cvd_reactor.set_species_floor(-1.0e-10)

#####################################################
# Run the inlet precursor composition parameter study
# ===================================================
# In the parameter study, the composition of the inlet SiF\ :sub:`4`\ molar ratio
# is done by resetting the molar composition of the inlet stream
# ``precursors`` with the recipes of the ``inlet_recipes`` list.
#
# Use the ``process_solution()`` method to convert the result from each PSR run
# into a solution ``Stream``. The surface site fractions and the bulk growth rates
# can be extract from the solution stream ``solnmixture`` with the
# ``surface_chemistry`` methods ``get_all_site_frac()`` and
# ``get_all_bulk_growth_rates()``, respectively.
# You can use the ``rop_surf()`` method to get the surface ROP of the gas
# species, the site species, and the bulk species of a surface material.
#
# .. note ::
#   In Pychemkin, the surface material is referred by its name. Instead of
#   the surface material index, you loop over the surface material names.
#   Use ``get_material_names()`` method to get a list of surface material
#   names defined in the surface mechanism.
#

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

##################################
# Plot the parameter study results
# ================================
# Plot the steady-state PSR solution against the inlet SiF\ :sub:`4` mole fraction.
# You should see that both the Si(D) growth rate and the NH\ :sub:`3` surface
# *consumption rate* increase with the increasing inlet SiF\ :sub:`4` mole fraction.
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
