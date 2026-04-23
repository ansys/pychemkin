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

r""".. _ref_trans_PSR_ignition:

===========================================================
Modeling atomic layer deposition process in a transient PSR
===========================================================

Atomic layer deposition (ALD) is a technique used to deposit thin films
of solid materials in a very controlled manner, and differs from CVD primarily
in that it is a transient process with the deposition surface being exposed to
pulses of alternating gases. Ideally, the deposition chemistry in ALD is
self-limiting, with growth occurring in a layer-by-layer manner and the deposition
thickness being controlled only by the number of cycles.

The major advantage of ALD over CVD is the improved control over the deposition
process and more conformal deposition. The inherently lower deposition rates,
however, lead to longer process times and higher costs. During a pulse, it is
important that enough molecules react with all parts of the substrate to be coated.
But many of the precursor materials are expensive. Thus one of the process
optimization goals is to reduce the amount of precursor that flows through the
reactor without reacting at the surface.

Normally the perfectly stirred reactor (PSR) model is a 0-D steady-state model with
constant reactor pressure. When the dynamic responses of a chemistry dominated process
are of interest, performing transient simulations of the system becomes necessary.
The transient PSR (transPSR) model of the Ansys Chemkin is a useful tool for the ALD
process development. It can be applied for optimizing the pulse sequences
and the cycle times of the ALD process to obtain the desired pattern and thickness of
the deposit layers. It can also be used to develop and to validate the chemical
reaction mechanism for the ALD process.

This project demonstrates the process of setting up and running the transPSR model
with multiple inlets of different flow rate profiles. The alumina ALD chemistry used
in this example is described in section 5.4.4 (Alumina ALD) of the Chemkin Tutorial
manual.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_trans_PSR_ALD.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger

# Chemkin 0-D transient PSR model with given reactor volume
from ansys.chemkin.core.stirreactors.transient_PSR import (
    TransientPSRSetVolumeFixedTemperature as TransPsr,
)
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

##########################################
# Create an instance of the Chemistry Set
# ========================================
# The chemistry set describes the atomic layer deposition (ALD) of
# alumina from trimethylaluminum (TMA) and ozone. This mechanism is
# deliberately simplistic for illustration purposes only. It demonstrates
# one way of describing the ALD of alumina, but it should generally be
# considered as illustrative only and not used as a source of kinetic data
# for this process. This mechanism is designed to deposit stoichiometric alumina,
# with three oxygen atoms being deposited for each two aluminum atoms.
#
# The ALD mechanism is provided in the *examples\data* folder of the local
# *PyChemkin* installation. You must provide the mechanism file names,
# ``chem_AL2O3_ald.inp`` and ``surf_AL2O3_ald.inp``, with the correct path
# ``example_mech_data`` to the Chemistry Set before pre-processing. The thermodynamic
# data file from the default Ansys Chemkin installation will be used.
#

# set mechanism directory (the default Chemkin mechanism data directory)
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
MyALDMech = ck.Chemistry(label="Al2O3_ALD")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "AL2O3_ALD"
MyALDMech.chemfile = str(example_mech_data / local_data_dir / "chem_AL2O3_ald.inp")
MyALDMech.surffile = str(example_mech_data / local_data_dir / "surf_AL2O3_ald.inp")
MyALDMech.thermfile = str(data_dir / "therm.dat")

##############################
# Preprocess the Chemistry Set
# ============================

# preprocess the mechanism files
ierror = MyALDMech.preprocess()
if ierror != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {ierror}")
    exit()

##########################
# Set up the inlet streams
# ========================
# Instantiate a stream feed for the inlet gas mixture.
# The stream is a mixture with the addition of the
# inlet flow rate. You specify inlet gas properties in the same way as you
# set up a mixture.
#
# Use the ``x`` method of the ``Stream`` object to specify the inlet gas
# composition with a ``recipe``. In this project, the ``mass_flowrate`` method,
# [g/sec] is used set the volumetric flow rates of the two inlet streams.
# The ``air`` composition is copied from the pre-defined ``Air`` object in Pychemkin.

# set the TMA stream composition and conditions
met_organic = Stream(MyALDMech, label="TMA")
met_organic.pressure = 1.0 * ck.P_TORRS  # 1 [torr]
met_organic.temperature = 620.2  # [K]
met_organic.x = [("ALMe3", 0.01), ("AR", 0.99)]
# TMA stream sccm flow rate profile in pulses
# TMA stream sccm data points
tma_time = np.array(
    [0.0,
     0.19, 
     0.2,
     5.19,
     5.2,
     5.39,
     5.4,
     10.39,
     10.4,
     10.59,
     10.6,
     15.59,
     15.6,
     15.79,
     15.8,
     20.8]
)  # [sec]
tma_profile = np.array(
    [800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0]
)  # [sccm]
# set the TMA stream sccm profile
met_organic.set_sccm_flowrate_profile(tma_time, tma_profile)

# set the air composition and conditions
oxidizers = Stream(MyALDMech, label="O3")
oxidizers.pressure = met_organic.pressure
oxidizers.temperature = 620.2
oxidizers.x = [("O2", 0.4), ("O3", 0.05), ("AR", 0.55)]
# oxidizers stream sccm flow rate profile in pulses
# oxidizers stream sccm data points
oxid_time = np.array(
    [0.0,
     1.19,
     1.2,
     3.19,
     3.2,
     6.39,
     6.4,
     8.39,
     8.4,
     11.59,
     11.6,
     13.59,
     13.6,
     16.79,
     16.8,
     18.79,
     18.8,
     20.8]
)  # [sec]
oxid_profile = np.array(
    [0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0]
)  # [sccm]
# set the oxidizers stream sccm profile
oxidizers.set_sccm_flowrate_profile(oxid_time, oxid_profile)

# set the air composition and conditions
purge = Stream(MyALDMech, label="AR")
purge.pressure = met_organic.pressure
purge.temperature = 620.2
purge.x = [("AR", 1.0)]
# purge stream sccm flow rate profile in pulses
# purge stream sccm data points
purge_time = np.array(
    [0.0,
     0.19,
     0.2,
     1.19,
     1.2,
     3.19,
     3.2,
     5.19,
     5.2,
     5.39,
     5.4,
     6.39,
     6.4,
     8.39,
     8.4,
     10.39,
     10.4,
     10.59,
     10.6,
     11.59,
     11.6,
     13.59,
     13.6,
     15.59,
     15.6,
     15.79,
     15.8,
     16.79,
     16.8,
     18.79,
     18.8,
     20.8]
)  # [sec]
purge_profile = np.array(
    [0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0,
     0.0,
     0.0,
     800.0,
     800.0]
)  # [sccm]
# set the purge stream sccm profile
purge.set_sccm_flowrate_profile(purge_time, purge_profile)
# check the initial flow rates of the streams
t = 0.0
print("Initial 'TMA' stream flow rate = "
      f"{met_organic.get_sccm_flowrate_profile_data(time=t)} [SCCM]")
print("Initial 'O3' stream flow rate = "
      f"{oxidizers.get_sccm_flowrate_profile_data(time=t)} [SCCM]")
print("Initial 'AR purge' stream flow rate = "
      f"{met_organic.get_sccm_flowrate_profile_data(time=t)} [SCCM]")

#######################################
# Set up the transient PSR bomb reactor
# =====================================
# Instantiate the transient PSR ``ald_reactor`` as a
# ``TransientPSRSetVolumeFixedTemperature`` object.
# Since this is a transient simulation, a initial reactor condition is required.
# In this example, the ``ald_reactor`` is initially filled with pure argon.

ald_reactor = TransPsr(purge, label="ALD_reactor")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For a PSR, you must provide
# either the residence time or the reactor volume. You can also make changes to
# any reactor initial conditions such as the reactor temperature and the gas
# composition. For example, the ``reset_initial_temperature()``` method and
# the ``reset_initial_composition()`` method.
# You can use the solver parameters to improve the convergence performance.

# set PSR volume (cm3): required for the
# ``TransientPSRSetVolumeFixedTemperature`` model
ald_reactor.volume = 600.0
# reactor active surface area [cm2]
ald_reactor.area = 300.0
# transient simulation end time [sec]
ald_reactor.time = 20.8
# initial surface conditions
# get the number of surface material defined in the surface mechanism
n_material = ald_reactor.get_numb_material
# get all surface material names
ald_material_names = ald_reactor.get_material_names()
# set up surface coverage of all surface materials
# initial surface coverage (site fractions)
site_recipe = [[("O(S)", 1.0)]]
# bulk activities
bulk_recipe = [[("AL2O3(B)", 1.0)]]
# set the surface condition of the surface material
for i, mname in enumerate(ald_material_names):
    # set area fraction [-] (optional: by default, materials have the same area)
    ald_reactor.set_material_area_fraction(mname, 1.0)
    # set material temperature [K] if it is different from the gas temperature
    # (optional: by default, it is the same as the mixture temperature)
    ald_reactor.set_material_temperature(mname, 723.0)
    # set reactor surface coverage of the material
    # site fractions (optional: by default, all sites have the same fraction)
    ald_reactor.set_site_fraction(mname, site_recipe[i])
    # bulk activities (optional: by default, all bulk species activities = 1.0)
    ald_reactor.set_bulk_activity(mname, bulk_recipe[i])

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

# connect the metal organic stream to the reactor
ald_reactor.set_inlet(met_organic)

# connect the oxidizer stream to the reactor
ald_reactor.set_inlet(oxidizers)

# connect the purge stream to the reactor
ald_reactor.set_inlet(purge)

# check the net inlet flow rate to the ALD reactor at time = 0.0 [sec]
# the net mass flow rate to the PSR must not be zero at any given time
# during the simulation
print(f"initial net mass flow rate to the reactor = "
      f"{ald_reactor.net_mass_flowrate} [g/sec]")

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set solver the tolerances
ald_reactor.tolerances = (1.0e-20, 1.0e-8)
# solver non-negative mode
ald_reactor.force_nonnegative = True
# transient solver max time step size [sec]
# ald_reactor.set_solver_max_timestep_size(2.0e-3)
# set adaptive solution saving distance
ald_reactor.adaptive_solution_saving(mode=True, steps=20)
# use the legacy solver to get better coonvergence performance for
# small mechanism
# ald_reactor.set_legacy_option(option=True)

#############################
# Run the transient PSR model
# ============================
# Once all the necessary input parameters are provided, use the ``run()`` method
# to start the transient PSR simulation. After the simulation
# is completed successfully, Use the ``process_solution()`` method to retrieve the raw
# solution profiles. If you are interested in the exit solution, use the
# ``get_last_solution_mixture()`` method to get the solution at the exit as a
# ``Stream`` object ``combustor_exhaust``.
#
# .. note ::
#   In Pychemkin, the surface material is referred by its name. Instead of
#   the surface material index, you loop over the surface material names.
#   Use ``get_material_names()`` method to get a list of surface material
#   names defined in the surface mechanism.
#
# .. note ::
#   For a site or a bulk species, use any ``Mixture`` or ``Stream`` created
#   with the ``Chemistry Set`` containing the surface mechanism to get its
#   species index. There are two type of index for surface species, the
#   global index includes all species (the gas plus the site and the bulk of
#   all materials); the local index includes the same type of surface species
#   (site or bulk) of all materials.
#


# set the start wall time
start_time = time.time()
# run the reactor model
runstatus = ald_reactor.run()
#
if runstatus != 0:
    # run failed
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)

runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec].")

# postprocess the solution profiles
ald_reactor.process_solution()

# get the number of data points in the solution profiles
solutionpoints = ald_reactor.getnumbersolutionpoints()
print(f"number of solution time points = {solutionpoints}")
# store solution profiles
# simulation time [sec]
sim_time = ald_reactor.get_solution_variable_profile("time")
# solution arrays
# get all surface material names
aldt_material_names = ald_reactor.get_material_names()
# AL2O3(B) linear growth rate [micron/min]
al2o3_b_growth_solution = np.zeros(solutionpoints, dtype=np.double)
al2o3_b_thickness = np.zeros_like(al2o3_b_growth_solution, dtype=np.double)
# get gas-species index
i_o = met_organic.get_specindex("O")
i_tma = met_organic.get_specindex("ALMe3")
o_solution = np.zeros_like(al2o3_b_growth_solution, dtype=np.double)
tma_solution = np.zeros_like(al2o3_b_growth_solution, dtype=np.double)
# get surface site fraction solution profiles
# get surface species index
# global index includes all species: the gas, the site of all materials, and the bulk
# of all materials
# local index includes the same type of species (site or bulk) of all materials
i_o_s_global, i_o_s_local = met_organic.get_surf_specindex("O(S)")
i_al2o3_b_global, i_al2o3_b_local = met_organic.get_surf_specindex("AL2O3(B)")
# surface O(S) site fraction
o_s_solution = ald_reactor.get_solution_variable_profile("O(S)")
# surface ALMe2(S) site fraction
alme2_s_solution = ald_reactor.get_solution_variable_profile("ALMe2(S)")
# surface ALMeOALMe(S) site fraction
almeoalme_s_solution = ald_reactor.get_solution_variable_profile("ALMeOALMe(S)")
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = ald_reactor.get_solution_mixture_at_index(solution_index=i)
    # calculate species production rate due to surface reactions [mole/cm2-sec] 
    brate = 0.0
    for m in ald_material_names:
        # calculate species ROPs due to surface reactions [mole/cm2-sec] 
        surf_rop, _ = solutionmixture.rop_surf(m)
        # get AL2O3(B) linear growth rate [cm/sec] due to surface reactions
        brates = (
            solutionmixture.surface_chemistry.get_bulk_linear_growth_rates(
                m, surf_rop
            )
        )
        brate += brates[i_al2o3_b_local]
    # linear growth rate of AL2O3(B) bulk species [cm/sec] -> [micron/sec]
    al2o3_b_growth_solution[i] = brate * 1.0e4
    # compute deposit thickness [micron]
    if i > 0:
        im = i - 1
        delta_time = sim_time[i] - sim_time[im]
        delta_grate = (al2o3_b_growth_solution[im] + al2o3_b_growth_solution[i])
        # integration
        al2o3_b_thickness[i] = (al2o3_b_thickness[im]
                                + delta_grate * delta_time / 2.0)
    # gas phase O solution profile [mole fraction]
    o_solution[i] = solutionmixture.x[i_o]
    # gas phase ALMe3 solution profile [mole fraction]
    tma_solution[i] = solutionmixture.x[i_tma]

# clean up
ck.done()

##########################################
# Plot the transient PSR solution profiles
# ========================================
# Plot the solution profiles over the entire simulation time span.
# You could observe how the growth of the AL2O3(B) thickness is controlled
# by turning on and off the inlet streams.

plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Transient PSR ALD Solution", fontsize=16)
plt.subplot(221)
plt.plot(tma_time, tma_profile, "b-")
plt.plot(oxid_time, oxid_profile, "r-")
plt.legend(["met_organic", "oxidizers"], loc="upper right")
plt.ylabel("Inlet volumetric flow rate [SCCM]")
plt.subplot(222)
plt.plot(sim_time, o_solution, "r-")
plt.ylabel("O Mole fraction")
plt.subplot(223)
plt.plot(sim_time, o_s_solution, "r-")
plt.plot(sim_time, almeoalme_s_solution, "b--")
plt.plot(sim_time, alme2_s_solution, "g:")
plt.legend(["O(S)", "ALMeOALMe(S)", "ALMe2(S)"], loc="upper right")
plt.xlabel("Time [sec]")
plt.ylabel("Site Fraction")
plt.subplot(224)
plt.plot(sim_time, al2o3_b_thickness, "b-")
plt.xlabel("Time [sec]")
plt.ylabel("Deposit thickness AL2O3(B) [micron]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_trans_PSR_ALD.png", bbox_inches="tight")
