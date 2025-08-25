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
.. _ref_rcm_multiproc:

====================================
Simulate a rapid compression machine
====================================

**This example runs multiple RCM simulation on several processors. The conditions are read from
a .txt file, including species concentrations, temperature, pressure, and time-volume profiles
The compression time is determined from the provided time-volume profile.
**The results, including the ignition delay times, temperature and pressure after compression are
collected, plotted, and saved.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_RCM_mul_solution.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
# chemkin batch reactor models (transient)
from ansys.chemkin.batchreactors.batchreactor import (
    GivenVolumeBatchReactor_EnergyConservation,
)
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching
import multiprocessing as mp  #multi processors

#################################
# initial conditions of RCM cases
# ===============================
# The initial conditions for RCM cases are provided in a txt file
# For nth case:
# Line 4n + 1: Initial species mole fractions
# Line 4n + 2: Initial temperature (K) and initial pressure (atm)
# Line 4n + 3 and 4n + 4: Time (s) and volume profile
# All numbers are separated by space, tab, comma, or their combinations (empty lines will be ignored)

# initial conditions for one single RCM case
class RCMcondition:
    def __init__(self, ind, speciesList, speciesFrac, temperature, pressure, time, volumeProfile):
        self.id = ind
        self.speciesList = speciesList
        self.speciesFrac = speciesFrac
        self.temperature = temperature
        self.pressure = pressure
        self.time = time
        self.volumeProfile = volumeProfile
        # determine the compression time from the volume profiles
        compressEnd = self.volumeProfile.index(min(self.volumeProfile))
        self.compressionTime = self.time[compressEnd]

    # this sets the composition, temperature, and pressure of a mixture
    def setRCMmixutre(self, premixed):
        if not isinstance(premixed, ck.Mixture):
            print(
                Color.RED + "** the first argument must be a Mixture object",
                end=Color.END,
            )
            raise TypeError
        premixed.X = list(zip(self.speciesList, self.speciesFrac))
        premixed.temperature = self.temperature
        premixed.pressure = self.pressure * ck.Patm

# initial conditions for multiple RCM cases,
class RCMconditions:
    def __init__(self, species, filename):
        self.file = filename
        self.species = species
        self.numSpecies = len(species)
        self.conditions = []
        i = 0
        # read conditions from file into list
        with open(self.file, "r") as file:
            for line in file:
                if line.strip():
                    cleaned_line = line.replace(',', ' ').replace('\t', ' ')
                    row = [float(num) for num in cleaned_line.split() if num.strip()]
                    if i % 4 == 0:
                        if len(row) != self.numSpecies:
                            print("Line #", i+1, "has wrong number of species")
                            raise ValueError
                        speciesFrac = row
                    elif i % 4 == 1:
                        if len(row) != 2:
                            print("Line #", i+1, "is expected to have 2 numbers")
                            raise ValueError
                        temperature = row[0]
                        pressure = row[1]
                    elif i % 4 == 2:
                        time = row
                    else:
                        volumeProfile = row
                        if len(row) != len(time):
                            print("Line #", i+1, "number of volume values does not match that of time values")
                            raise ValueError
                        condition = RCMcondition(i//4+1, species, speciesFrac, temperature, pressure, time, volumeProfile)
                        self.conditions.append(condition)
                    i = i + 1
        if i % 4 != 0:
            print("Conditions incomplete, missing lines")
            raise ValueError
        self.numConditions = i // 4

###################
#  Run one RCM case
# =================
# This function runs one RCM case at given path with the mechanism and initial conditions provided

def run_rcm(work_path, chemfile, thermfile, condition):
    # set current work directory
    os.chdir(work_path)
    # set verbose mode
    ck.set_verbose(True)

#######################################
# Create and preprocess a chemistry set
# =====================================

    MyGasMech = ck.Chemistry(label="GRI 3.0")
    # set mechanism input files
    MyGasMech.chemfile = chemfile
    MyGasMech.thermfile = thermfile
    # preprocess the mechanism files
    iError = MyGasMech.preprocess()

######################################
# Set up the rapid-compression machine
# ====================================

    # create the premixed mixture
    premixed = ck.Mixture(MyGasMech)
    # pass initial conditions to the mixture
    condition.setRCMmixutre(premixed)
    # create a constant volume batch reactor (with energy equation)
    MyCONV = GivenVolumeBatchReactor_EnergyConservation(premixed, label="RCM")
    # set volume profile
    MyCONV.set_volume_profile(condition.time, condition.volumeProfile)

############################################
# Set up additional reactor model parameters
# ==========================================
# You must provide reactor parameters, solver controls, and output instructions
# before running the simulations. For a batch reactor, the initial volume and the
# simulation end time are required inputs.

    MyCONV.time = 0.8  # sec

####################
# Set output options
# ==================

    # set timestep between saving solution
    MyCONV.timestepforsavingsolution = 0.001
    # turn ON saving to XML solution file (default)
    MyCONV.XML_Output = True
    # turn ON adaptive solution saving
    MyCONV.adaptive_solution_saving(mode=True, value_change=100, target="TEMPERATURE")
    # change timestep between saving solution
    MyCONV.timestepforsavingsolution = 0.01
    # set ignition delay criteria
    MyCONV.set_ignition_delay("T_rise", 600)
    #MyCONV.set_ignition_delay(method="Species_peak", target="OH")

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods, such as those
# for tolerances.

    # set tolerance
    MyCONV.tolerances = (1.0e-10, 1.0e-8)
    # turn on the force non-negative solutions option in the solver
    MyCONV.forcenonnegative = True

####################
# run the simulation
#===================
# Use the ``run()`` method to start the RCM simulation.

    # run the CONV reactor model
    runstatus = MyCONV.run()
    # check run status
    if runstatus != 0:
        # run failed!
        print(Color.RED + f">>> Case {condition.id} RUN FAILED <<<", end=Color.END)
        exit()
    # run success!
    print(Color.GREEN + f">>> Case {condition.id} RUN COMPLETED <<<", end=Color.END)

######################################
# get and plot the temperature profile
#=====================================

    # get ignition delay time and change the unit to ms
    ignitionDelaytime = round(MyCONV.get_ignition_delay() - condition.compressionTime * 1000, 2)
    # post-process the solutions
    MyCONV.process_solution()
    # get time-temperature profile
    timeprofile = MyCONV.get_solution_variable_profile("time")
    tempprofile = MyCONV.get_solution_variable_profile("temperature")
    # get temperature and pressure at the end of compression
    compressTemp = MyCONV.get_solution_mixture(condition.compressionTime).temperature
    compressPres = MyCONV.get_solution_mixture(condition.compressionTime).pressure / ck.Patm
    # plot and save time-temperature profile
    parent_dir = os.path.abspath(os.path.join(work_path, os.pardir))
    plt.plot(timeprofile, tempprofile, "r-")
    plt.ylabel("Temperature [K]")
    plt.savefig(os.path.join(parent_dir, f"TemperatureProfile_{condition.id}.png"), format="png", dpi=200)
    plt.clf()
    return ignitionDelaytime, compressTemp, compressPres

#########################
#  Run multiple RCM cases
# =======================
# This function runs multiple RCM cases over a specified number of processors
# initial conditions are provided in a txt file

def run_rcm_multiproc():

    # The mechanism to load is the GRI 3.0 mechanism for methane combustion.
    # This mechanism and its associated data files come with the standard Ansys Chemkin
    # installation in the ``/reaction/data`` directory.
    mechanism_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
    thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

    # check working directory
    current_dir = os.getcwd()
    print("current dir is" + current_dir)
    # set number of processors
    numProcessors = 4
    # specify txt file name for initial conditions
    conditionFile = "data/rcm_multiproc_conditions.txt"
    # read RCM conditions
    rcmInput = RCMconditions(["CH4", "N2", "O2"], os.path.join(current_dir, conditionFile))

#############################################
# run the simulation with multiple processors
# ===========================================
# Use ``mp.Pool()``, ``run_rcm()`` function and ``args`` to run simulations with multiple processors

    args = []
    # make a new directory for each case and run all cases over multiple processors
    for i in range(rcmInput.numConditions):
        work_path = os.path.join(current_dir, "RUN_" + str(i+1).zfill(4))
        if not os.path.exists(work_path):
            os.makedirs(work_path)
        args.append((work_path, chemfile, thermfile, rcmInput.conditions[i]))
    with mp.Pool(processes=numProcessors) as pool:
        results = pool.starmap(run_rcm, args)

#################################
# collect, plot, and save results
#================================

    # collect results
    ignitionDelaytimes = []
    compressTemp = []
    compressPres = []
    for i in range(rcmInput.numConditions):
        ignitionDelaytimes.append(results[i][0])
        compressTemp.append(results[i][1])
        compressPres.append(results[i][2])
    # plot results
    reciprocalT = 1000/np.array(compressTemp)
    plt.plot(reciprocalT, ignitionDelaytimes, "b-")
    plt.yscale("log")
    plt.xlabel("Time [s]")
    plt.ylabel("Ignition delay times [ms]")
    plt.xlabel("1000/T [K]")
    plt.savefig("plot_RCM_mul_solution.png", format="png", dpi=200)
    plt.show()
    # save results to a txt file
    ids = np.arange(1, rcmInput.numConditions + 1)
    resultsRCM = np.column_stack((ids, compressPres, compressTemp, ignitionDelaytimes))
    np.savetxt("ignition_delaytimes.txt", resultsRCM, fmt="%.1f", delimiter=" ", header="Pc(atm) Tc(K) IgnitionDelay(ms)")
    # print ignition delay times
    print(f"ignition delay times = {ignitionDelaytimes} [msec]")

if __name__ == "__main__":
    run_rcm_multiproc()
